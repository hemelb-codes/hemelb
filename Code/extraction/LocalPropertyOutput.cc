// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include "extraction/LocalPropertyOutput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/formats/offset.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "io/writers/xdr/XdrVectorWriter.h"
#include "net/IOCommunicator.h"
#include "constants.h"
#include "units.h"

namespace hemelb
{
  namespace extraction
  {
    namespace
    {
      // Declare recursive helper
      template <typename... Ts>
      io::writers::Writer& encode(io::writers::Writer& enc, Ts... args);
      // Terminating case - one arg
      template <typename T>
      io::writers::Writer& encode(io::writers::Writer& enc, T arg) {
	return enc << arg;
      }
      // Recursive case - N + 1 args
      template <typename T, typename... Ts>
      io::writers::Writer& encode(io::writers::Writer& enc, T arg, Ts... args) {
	return encode(enc << arg, args...);
      }

      // XDR encode some values and return the result buffer
      template <typename... Ts>
      std::vector<char> quick_encode(Ts... args) {
	io::writers::xdr::XdrVectorWriter encoder;
	encode(encoder, args...);
	auto ans = encoder.GetBuf();
	return ans;
      }

      // Helper for writing values converted to the type contained in
      // the code::Type variant tag value.
      //
      // Use of std::variant + visit ensures that we generate all the
      // types required with only a single implementation.
      template <typename XDRW, typename... MemTs>
      void write(XDRW& writer, code::Type tc, MemTs... vals) {
	using common_t = std::common_type_t<MemTs...>;
	static_assert(std::conjunction_v<std::is_same<common_t, MemTs>...>,
		      "Values must be of uniform type");

	std::visit([&] (auto tag) {
	    using FileT = decltype(tag);
	    (writer << ... << FileT(vals));
	  },
	  tc);
      }

    }  // namespace

    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource,
                                             const PropertyOutputFile& outputSpec,
                                             const net::IOCommunicator& ioComms) :
        comms(ioComms), dataSource(dataSource), outputSpec(outputSpec)
    {
      // Open the file as write-only, create it if it doesn't exist, don't create if the file
      // already exists.
      outputFile = net::MpiFile::Open(comms, outputSpec.filename,
                                      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);

      // Create a new file name by changing the suffix - outputSpec
      // must have a .-separated extension!
      const auto offsetFileName = io::formats::offset::ExtractionToOffset(outputSpec.filename);

      // Now create the file.
      offsetFile = net::MpiFile::Open(comms, offsetFileName,
				      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);

      // Count sites on this task
      uint64_t siteCount = 0;
      dataSource.Reset();
      while (dataSource.ReadNext())
      {
        if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition()))
        {
          ++siteCount;
        }
      }

      // Calculate how long local writes need to be.

      // First get the length per-site
      // Always have 3 uint32's for the position of a site
      writeLength = 3 * 4;

      // Then get add each field's length
      for (auto&& f: outputSpec.fields) {
	// Also check that len offsets makes sense
	auto n = f.noffsets;
	auto len = GetFieldLength(f.src);
	if (n == 0 || n == 1 || n == len) {
	  // ok
	} else {
	  throw Exception() << "Invalid length of offsets array " << n;
	}
        writeLength += len * code::type_to_size(f.typecode);
      }

      //  Now multiply by local site count
      writeLength *= siteCount;

      // The IO proc also writes the iteration number
      if (comms.OnIORank())
      {
        writeLength += 8;
      }

      //! @TODO: These two MPI calls can be replaced with one

      // Everyone needs to know the total length written during one iteration.
      allCoresWriteLength = comms.AllReduce(writeLength, MPI_SUM);

      // Only the root process needs to know the total number of sites written
      // Note this has a garbage value on other procs.
      uint64_t allSiteCount = comms.Reduce(siteCount, MPI_SUM, comms.GetIORank());

      unsigned totalHeaderLength = 0;

      // Write the header information on the IO proc.
      if (comms.OnIORank())
      {
        // Compute the length of the field header
        unsigned const fieldHeaderLength = std::transform_reduce(
	  outputSpec.fields.begin(), outputSpec.fields.end(),
	  0U,
	  std::plus<unsigned>{},
	  [&](OutputField const& f) {
	    return io::formats::extraction::GetFieldHeaderLength(
	      f.name,
	      f.noffsets,
	      code::type_to_enum(f.typecode)
	    );
	  }
	);

        totalHeaderLength = io::formats::extraction::MainHeaderLength + fieldHeaderLength;

	io::writers::xdr::XdrVectorWriter headerWriter;

	// Encoder for ONLY the main header (note shorter length)
	headerWriter << uint32_t(io::formats::HemeLbMagicNumber)
		     << uint32_t(io::formats::extraction::MagicNumber)
		     << uint32_t(io::formats::extraction::VersionNumber);
	headerWriter << double(dataSource.GetVoxelSize());
	const util::Vector3D<distribn_t> &origin = dataSource.GetOrigin();
	headerWriter << double(origin[0]) << double(origin[1]) << double(origin[2]);

	// Write the total site count and number of fields
	headerWriter << uint64_t(allSiteCount) << uint32_t(outputSpec.fields.size())
		     << uint32_t(fieldHeaderLength);

	// Main header now finished - do field headers
	for (auto& field: outputSpec.fields) {
	  auto const len = GetFieldLength(field.src);
	  headerWriter << field.name
		       << uint32_t(len)
		       << uint32_t(code::type_to_enum(field.typecode))
		       << field.noffsets;
	  std::visit([&](auto&& tag) {
	      for(auto& offset: field.offset)
		headerWriter << (decltype(tag))offset;
	    },
	    field.typecode);
	}

	assert(headerWriter.GetBuf().size() == totalHeaderLength);
        // Write from the buffer
        outputFile.WriteAt(0, headerWriter.GetBuf());
      }

      // Calculate where each core should start writing
      if (comms.OnIORank())
      {
        // For core 0 this is easy: it passes the value for core 1 to the core.
        localDataOffsetIntoFile = totalHeaderLength;

        if (comms.Size() > 1)
        {
          comms.Send(localDataOffsetIntoFile + writeLength, 1, 1);
        }
      }
      else
      {
        // Receive the writing start position from the previous core.
        comms.Receive(localDataOffsetIntoFile, comms.Rank() - 1, 1);

        // Send the next core its start position.
        if (comms.Rank() != (comms.Size() - 1))
        {
          comms.Send(localDataOffsetIntoFile + writeLength, comms.Rank() + 1, 1);
        }
      }

      // Create the buffer that we'll write each iteration's data into.
      buffer.resize(writeLength);

      WriteOffsetFile();
    }

    LocalPropertyOutput::~LocalPropertyOutput()
    {

    }

    bool LocalPropertyOutput::ShouldWrite(unsigned long timestepNumber) const
    {
      return ( (timestepNumber % outputSpec.frequency) == 0);
    }

    const PropertyOutputFile& LocalPropertyOutput::GetOutputSpec() const
    {
      return outputSpec;
    }

    void LocalPropertyOutput::Write(unsigned long timestepNumber)
    {
      // Don't write if we shouldn't this iteration.
      if (!ShouldWrite(timestepNumber))
      {
        return;
      }

      // Don't write if this core doesn't do anything.
      if (writeLength <= 0)
      {
        return;
      }

      // Create the buffer.
      auto xdrWriter = io::MakeXdrWriter(buffer.begin(), buffer.end());

      // Firstly, the IO proc must write the iteration number.
      if (comms.OnIORank())
      {
        xdrWriter << (uint64_t) timestepNumber;
      }

      dataSource.Reset();

      while (dataSource.ReadNext())
      {
        const util::Vector3D<site_t>& position = dataSource.GetPosition();
        if (outputSpec.geometry->Include(dataSource, position))
        {
          // Write the position
          xdrWriter << (uint32_t) position.x << (uint32_t) position.y << (uint32_t) position.z;

          // Write for each field.
          for (auto& fieldSpec: outputSpec.fields)
          {
	    overload_visit(
	      fieldSpec.src,
	      [&](source::Pressure) {
		write(xdrWriter, fieldSpec.typecode, dataSource.GetPressure() - fieldSpec.offset[0]);
	      },
              [&](source::Velocity) {
		auto&& v = dataSource.GetVelocity();
		write(xdrWriter, fieldSpec.typecode, v.x, v.y, v.z);
	      },
	      //! @TODO: Work out how to handle the different stresses.
              [&](source::VonMisesStress) {
		write(xdrWriter, fieldSpec.typecode, dataSource.GetVonMisesStress());
              },
              [&](source::ShearStress) {
		write(xdrWriter, fieldSpec.typecode, dataSource.GetShearStress());
	      },
              [&](source::ShearRate) {
		write(xdrWriter, fieldSpec.typecode, dataSource.GetShearRate());
	      },
              [&](source::StressTensor) {
		util::Matrix3D tensor = dataSource.GetStressTensor();
		// Only the upper triangular part of the symmetric tensor is stored. Storage is row-wise.
                write(xdrWriter, fieldSpec.typecode,
                      tensor[0][0], tensor[0][1], tensor[0][2],
                                    tensor[1][1], tensor[1][2],
                                                  tensor[2][2]);
              },
              [&](source::Traction) {
		auto&& t = dataSource.GetTraction();
		write(xdrWriter, fieldSpec.typecode, t.x, t.y, t.z);
	      },
              [&](source::TangentialProjectionTraction) {
		auto&& t = dataSource.GetTangentialProjectionTraction();
		write(xdrWriter, fieldSpec.typecode, t.x, t.y, t.z);
	      },
              [&](source::Distributions) {
		unsigned numComponents= dataSource.GetNumVectors();
		distribn_t const* d_ptr = dataSource.GetDistribution();
		for (int i = 0; i < numComponents; i++)
                {
		  write(xdrWriter, fieldSpec.typecode, d_ptr[i]);
		}
	      },
              [&](source::MpiRank) {
                write(xdrWriter, fieldSpec.typecode, comms.Rank());
	      }
            );
          }
        }
      }

      // Actually do the MPI writing.
      outputFile.WriteAt(localDataOffsetIntoFile, buffer);

      // Set the offset to the right place for writing on the next iteration.
      localDataOffsetIntoFile += allCoresWriteLength;
    }

    // Write the offset file.
    void LocalPropertyOutput::WriteOffsetFile() {
      namespace fmt = io::formats;

      // On process 0 only, write the header
      if (comms.OnIORank()) {
	auto buf = quick_encode(
				uint32_t(fmt::HemeLbMagicNumber),
				uint32_t(fmt::offset::MagicNumber),
				uint32_t(fmt::offset::VersionNumber),
				uint32_t(comms.Size())
				);
	assert(buf.size() == fmt::offset::HeaderLength);
	offsetFile.WriteAt(0, buf);
      }
      // Every rank writes its offset
      uint64_t offsetForOffset = comms.Rank() * sizeof(localDataOffsetIntoFile)
	+ fmt::offset::HeaderLength;
      offsetFile.WriteAt(offsetForOffset, quick_encode(localDataOffsetIntoFile));

      // Last process writes total
      if (comms.Rank() == (comms.Size()-1)) {
	offsetFile.WriteAt(offsetForOffset + sizeof(localDataOffsetIntoFile),
			   quick_encode(localDataOffsetIntoFile + writeLength));
      }
    }

    unsigned LocalPropertyOutput::GetFieldLength(source::Type src) const
    {
      return overload_visit(src,
	[](source::Pressure) {
	  return 1U;
	},
	[](source::Velocity) {
	  return 3U;
	},
	[](source::ShearStress) {
	  return 1U;
	},
	[](source::VonMisesStress) {
	  return 1U;
	},
	[](source::ShearRate) {
	  return 1U;
	},
	[](source::StressTensor) {
	  return 6U;
	},
	[](source::Traction) {
	  return 3U;
	},
	[](source::TangentialProjectionTraction) {
	  return 3U;
	},
	[&](source::Distributions) {
	  return dataSource.GetNumVectors();
	},
	[](source::MpiRank) {
	  return 1U;
	}
      );
    }
  }
}
