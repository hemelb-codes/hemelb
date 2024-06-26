// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "hassert.h"
#include "extraction/LocalPropertyOutput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/formats/offset.h"
#include "io/writers/XdrMemWriter.h"
#include "io/writers/XdrVectorWriter.h"
#include "net/IOCommunicator.h"
#include "util/span.h"
#include "constants.h"
#include "units.h"

namespace hemelb::extraction
{
    namespace
    {
      // Declare recursive helper
      template <typename... Ts>
      io::Writer& encode(io::Writer& enc, Ts... args);
      // Terminating case - one arg
      template <typename T>
      io::Writer& encode(io::Writer& enc, T arg) {
	return enc << arg;
      }
      // Recursive case - N + 1 args
      template <typename T, typename... Ts>
      io::Writer& encode(io::Writer& enc, T arg, Ts... args) {
	return encode(enc << arg, args...);
      }

      // XDR encode some values and return the result buffer
      template <typename... Ts>
      std::vector<char> quick_encode(Ts... args) {
	io::XdrVectorWriter encoder;
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

    static unsigned CalcFieldHeaderLength(std::vector<OutputField> const& fields);

    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource,
                                             const PropertyOutputFile& outputSpec_,
                                             const net::IOCommunicator& ioComms) :
        comms(ioComms), dataSource(dataSource), outputSpec(outputSpec_)
    {
      if (std::holds_alternative<multi_timestep_file>(outputSpec.ts_mode)) {
	// Just replace extension with .off
	offset_file_name = io::formats::offset::ExtractionToOffset(outputSpec.filename);
	// empty output_file_pattern is OK
      } else if (std::holds_alternative<single_timestep_files>(outputSpec.ts_mode)) {
	// Get views of the whole path
	std::string_view p = outputSpec.filename.native();
	auto i_pcd = p.find("%d", 0, 2);
	// The part before %d
	auto beginning = p.substr(0, i_pcd);
	// The part after
	auto end = p.substr(i_pcd + 2);
	// Construct the path without '%d'
	std::string basename{beginning};
	basename += end;
	// Use this to compute offset file name
	offset_file_name = io::formats::offset::ExtractionToOffset(basename);
	// Build the pattern
	output_file_pattern += beginning;
	output_file_pattern += "%*ld";
	output_file_pattern += end;
      }

      header_length = io::formats::extraction::MainHeaderLength + CalcFieldHeaderLength(outputSpec.fields);

      // Count sites on this rank
      local_site_count = CountWrittenSitesOnRank();
      global_site_count = comms.AllReduce(local_site_count, MPI_SUM);

      // Calculate how long local writes need to be (recall only IO
      // rank writes the timestep).
      auto const site_len = CalcSiteWriteLen(outputSpec.fields);
      local_data_write_length = local_site_count * site_len  + (comms.OnIORank() ? 8U : 0U);
      // Everyone needs to know the total length written during one iteration
      global_data_write_length = site_len * global_site_count + 8U;

      // Work out the offset for where this rank writes its data
      auto const local_write_end = comms.Scan(local_data_write_length, MPI_SUM) + header_length;
      local_write_start = local_write_end - local_data_write_length;

      // Prepare the header information on the IO proc.
      if (comms.OnIORank())
      {
	header_data = PrepareHeader();
      }

      // Create the buffer that we'll write each iteration's data into.
      buffer.resize(local_data_write_length);

      // Write the offset file
      WriteOffsetFile();

      // If we are doing all timesteps in one file, set it up now.
      if (std::holds_alternative<multi_timestep_file>(outputSpec.ts_mode)) {
	StartFile(outputSpec.filename);
      }
    }

    uint64_t LocalPropertyOutput::CountWrittenSitesOnRank() {
      auto n = uint64_t{0};
      dataSource.Reset();
      while (dataSource.ReadNext())
      {
	if (outputSpec.geometry->Include(dataSource, dataSource.GetPosition()))
        {
	  ++n;
	}
      }
      return n;
    }

    // Work out how many bytes are needed to write one site's data.
    std::uint64_t LocalPropertyOutput::CalcSiteWriteLen(std::vector<OutputField> const& fields) const {
      // Always have 3 uint32's for the position of a site
      std::uint64_t site_len = 3 * 4;

      // Then get add each field's length
      for (auto&& f: fields) {
	// Also check that len offsets makes sense
	auto n = f.noffsets;
	auto len = GetFieldLength(f.src);
	if (n == 0 || n == 1 || n == len) {
	  // ok
	} else {
	  throw Exception() << "Invalid length of offsets array " << n;
	}
        site_len += len * code::type_to_size(f.typecode);
      }
      return site_len;
    }

    // Compute the length of the field header
    unsigned CalcFieldHeaderLength(std::vector<OutputField> const& fields) {
      return std::transform_reduce(
        fields.begin(),fields.end(),
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
    }

    std::vector<char> LocalPropertyOutput::PrepareHeader() const {

      unsigned const field_header_len = CalcFieldHeaderLength(outputSpec.fields);
      unsigned const total_header_len = io::formats::extraction::MainHeaderLength + field_header_len;
      io::XdrVectorWriter headerWriter;

      // Encoder for ONLY the main header (note shorter length)
      headerWriter << std::uint32_t(io::formats::HemeLbMagicNumber)
		   << std::uint32_t(io::formats::extraction::MagicNumber)
		   << std::uint32_t(io::formats::extraction::VersionNumber);
      headerWriter << double(dataSource.GetVoxelSize());
      const util::Vector3D<distribn_t> &origin = dataSource.GetOrigin();
      headerWriter << double(origin[0]) << double(origin[1]) << double(origin[2]);

      // Write the total site count and number of fields
      headerWriter << std::uint64_t(global_site_count) << std::uint32_t(outputSpec.fields.size())
		   << std::uint32_t(field_header_len);

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

      HASSERT(headerWriter.GetBuf().size() == total_header_len);
      return headerWriter.GetBuf();
    }

    bool LocalPropertyOutput::ShouldWrite(unsigned long timestepNumber) const
    {
      return ( (timestepNumber % outputSpec.frequency) == 0);
    }

    const PropertyOutputFile& LocalPropertyOutput::GetOutputSpec() const
    {
      return outputSpec;
    }

    void LocalPropertyOutput::StartFile(std::string const& fn)
    {
      // Open the file as write-only, create it if it doesn't exist,
      // don't create if the file already exists.
      outputFile = net::MpiFile::Open(comms, fn,
                                      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);

      // Write the header information on the IO proc.
      if (comms.OnIORank())
      {
        // Write from the buffer
        outputFile.WriteAt(0, to_const_span(header_data));
      }
    }

    template <typename... Ts>
    std::string safe_fmt(std::string const& pattern, Ts... args) {
        int sz = std::snprintf(nullptr, 0,
                               pattern.data(), args...);
        if (sz < 0)
            throw Exception() << "Formatting error";

        // +1 for the null terminator
        std::string ans(sz + 1, '\0');
        std::snprintf(ans.data(), ans.size(),
                      pattern.data(), args...);
        return ans;
    }

    void LocalPropertyOutput::Write(unsigned long timestepNumber, unsigned long totalSteps)
    {
        // Don't write if we shouldn't this iteration.
        if (!ShouldWrite(timestepNumber))
        {
            return;
        }

        if (std::holds_alternative<single_timestep_files>(outputSpec.ts_mode)) {
            int prec = 3;
            unsigned long next = 1000;
            while (totalSteps > next) {
                prec += 1;
                next *= 10;
            }
            std::string fn = safe_fmt(output_file_pattern, prec, timestepNumber);
            StartFile(fn);
        }

      // Don't write if this core doesn't do anything.
      if (local_data_write_length > 0)
      {
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
	    xdrWriter << (uint32_t) position.x() << (uint32_t) position.y() << (uint32_t) position.z();

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
		  write(xdrWriter, fieldSpec.typecode, v.x(), v.y(), v.z());
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
		  // Only the upper triangular part of the symmetric
		  // tensor is stored. Storage is row-wise.
		  write(xdrWriter, fieldSpec.typecode,
			tensor[0][0], tensor[0][1], tensor[0][2],
                                      tensor[1][1], tensor[1][2],
                                                    tensor[2][2]);
		},
		[&](source::Traction) {
		  auto&& t = dataSource.GetTraction();
		  write(xdrWriter, fieldSpec.typecode, t.x(), t.y(), t.z());
		},
		[&](source::TangentialProjectionTraction) {
		  auto&& t = dataSource.GetTangentialProjectionTraction();
		  write(xdrWriter, fieldSpec.typecode, t.x(), t.y(), t.z());
		},
		[&](source::Distributions) {
		  unsigned numComponents = dataSource.GetNumVectors();
		  distribn_t const* d_ptr = dataSource.GetDistribution();
		  for (auto i = 0U; i < numComponents; i++)
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
	outputFile.WriteAt(local_write_start, to_const_span(buffer));
      }

      overload_visit(
        outputSpec.ts_mode,
	[this](multi_timestep_file) {
	  // Set the offset to the right place for writing on the next
	  // iteration.
	  local_write_start += global_data_write_length;
	},
	[this](single_timestep_files) {
	  outputFile.Close();
	}
      );
    }

    // Write the offset file.
    void LocalPropertyOutput::WriteOffsetFile() {
      namespace fmt = io::formats;

      // Create the file.
      auto offsetFile = net::MpiFile::Open(comms, offset_file_name,
				      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);

      // On process 0 only, write the header
      if (comms.OnIORank()) {
	auto buf = quick_encode(
				uint32_t(fmt::HemeLbMagicNumber),
				uint32_t(fmt::offset::MagicNumber),
				uint32_t(fmt::offset::VersionNumber),
				int32_t(comms.Size())
				);
	HASSERT(buf.size() == fmt::offset::HeaderLength);
	offsetFile.WriteAt(0, to_const_span(buf));
      }
      // Every rank writes its offset
      uint64_t offsetForOffset = comms.Rank() * sizeof(local_write_start)
	+ fmt::offset::HeaderLength;
      offsetFile.WriteAt(offsetForOffset, to_const_span(quick_encode(local_write_start)));

      // Last process writes total
      if (comms.Rank() == (comms.Size()-1)) {
	offsetFile.WriteAt(offsetForOffset + sizeof(local_write_start),
			   to_const_span(quick_encode(local_write_start + local_data_write_length)));
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
