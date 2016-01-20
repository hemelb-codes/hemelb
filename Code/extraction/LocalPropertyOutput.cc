
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include "extraction/LocalPropertyOutput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "net/IOCommunicator.h"
#include "constants.h"

namespace hemelb
{
  namespace extraction
  {
    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource,
                                             const PropertyOutputFile* outputSpec,
                                             const net::IOCommunicator& ioComms) :
      comms(ioComms), dataSource(dataSource), outputSpec(outputSpec)
    {
      // Open the file as write-only, create it if it doesn't exist, don't create if the file
      // already exists.
      outputFile = net::MpiFile::Open(comms, outputSpec->filename,
                                      MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);
      // Count sites on this task
      uint64_t siteCount = 0;
      dataSource.Reset();
      while (dataSource.ReadNext())
      {
        if (outputSpec->geometry->Include(dataSource, dataSource.GetPosition()))
        {
          ++siteCount;
        }
      }

      // Calculate how long local writes need to be.

      // First get the length per-site
      // Always have 3 uint32's for the position of a site
      writeLength = 3 * 4;

      // Then get add each field's length
      for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
      {
        writeLength += sizeof(WrittenDataType)
            * GetFieldLength(outputSpec->fields[outputNumber].type);
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
      uint64_t allSiteCount = comms.Reduce(siteCount, MPI_SUM,
                                           comms.GetIORank());

      unsigned totalHeaderLength = 0;

      // Write the header information on the IO proc.
      if (comms.OnIORank())
      {
        // Compute the length of the field header
        unsigned fieldHeaderLength = 0;
        for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
        {
          // Name
          fieldHeaderLength
              += io::formats::extraction::GetStoredLengthOfString(outputSpec->fields[outputNumber].name);
          // Uint32 for number of fields
          fieldHeaderLength += 4;
          // Double for the offset in each field
          fieldHeaderLength += 8;
        }

        // Create a header buffer
        totalHeaderLength = io::formats::extraction::MainHeaderLength + fieldHeaderLength;
        std::vector<char> headerBuffer(totalHeaderLength);

        {
          // Encoder for ONLY the main header (note shorter length)
          io::writers::xdr::XdrMemWriter
              mainHeaderWriter(&headerBuffer[0], io::formats::extraction::MainHeaderLength);

          // Fill it
          mainHeaderWriter << uint32_t(io::formats::HemeLbMagicNumber)
              << uint32_t(io::formats::extraction::MagicNumber)
              << uint32_t(io::formats::extraction::VersionNumber);
          mainHeaderWriter << double(dataSource.GetVoxelSize());
          const util::Vector3D<distribn_t> &origin = dataSource.GetOrigin();
          mainHeaderWriter << double(origin[0]) << double(origin[1]) << double(origin[2]);

          // Write the total site count and number of fields
          mainHeaderWriter << uint64_t(allSiteCount) << uint32_t(outputSpec->fields.size())
              << uint32_t(fieldHeaderLength);
          // Main header now finished.
          // Exiting the block kills the mainHeaderWriter.
        }
        {
          // Create the field header writer
          io::writers::xdr::XdrMemWriter
              fieldHeaderWriter(&headerBuffer[io::formats::extraction::MainHeaderLength],
                                fieldHeaderLength);
          // Write it
          for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
          {
            fieldHeaderWriter << outputSpec->fields[outputNumber].name
                << uint32_t(GetFieldLength(outputSpec->fields[outputNumber].type))
                << GetOffset(outputSpec->fields[outputNumber].type);
          }
          //Exiting the block cleans up the writer
        }

        // Write from the buffer
        outputFile.WriteAt(0, headerBuffer);
      }

      // Calculate where each core should start writing
      if (comms.OnIORank())
      {
        // For core 0 this is easy: it passes the value for core 1 to the core.
        localDataOffsetIntoFile = totalHeaderLength;

        if (comms.Size() > 1)
        {
          comms.Send(localDataOffsetIntoFile+writeLength, 1, 1);
        }
      }
      else
      {
        // Receive the writing start position from the previous core.
        comms.Receive(localDataOffsetIntoFile, comms.Rank()-1, 1);

        // Send the next core its start position.
        if (comms.Rank() != (comms.Size() - 1))
        {
          comms.Send(localDataOffsetIntoFile+writeLength, comms.Rank() + 1, 1);
        }
      }

      // Create the buffer that we'll write each iteration's data into.
      buffer.resize(writeLength);
    }

    LocalPropertyOutput::~LocalPropertyOutput()
    {

    }

    bool LocalPropertyOutput::ShouldWrite(unsigned long timestepNumber) const
    {
      return ( (timestepNumber % outputSpec->frequency) == 0);
    }

    const PropertyOutputFile* LocalPropertyOutput::GetOutputSpec() const
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
      io::writers::xdr::XdrMemWriter xdrWriter(&buffer[0], buffer.size());

      // Firstly, the IO proc must write the iteration number.
      if (comms.OnIORank())
      {
        xdrWriter << (uint64_t) timestepNumber;
      }

      dataSource.Reset();

      while (dataSource.ReadNext())
      {
        const util::Vector3D<site_t>& position = dataSource.GetPosition();
        if (outputSpec->geometry->Include(dataSource, position))
        {
          // Write the position
          xdrWriter << (uint32_t) position.x << (uint32_t) position.y << (uint32_t) position.z;

          // Write for each field.
          for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
          {
            switch (outputSpec->fields[outputNumber].type)
            {
              case OutputField::Pressure:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetPressure()
                    - REFERENCE_PRESSURE_mmHg);
                break;
              case OutputField::Velocity:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetVelocity().x)
                    << static_cast<WrittenDataType> (dataSource.GetVelocity().y)
                    << static_cast<WrittenDataType> (dataSource.GetVelocity().z);
                break;
                //! @TODO: Work out how to handle the different stresses.
              case OutputField::VonMisesStress:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetVonMisesStress());
                break;
              case OutputField::ShearStress:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetShearStress());
                break;
              case OutputField::ShearRate:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetShearRate());
                break;
              case OutputField::StressTensor:
              {
                util::Matrix3D tensor = dataSource.GetStressTensor();
                // Only the upper triangular part of the symmetric tensor is stored. Storage is row-wise.
                xdrWriter << static_cast<WrittenDataType> (tensor[0][0])
                    << static_cast<WrittenDataType> (tensor[0][1])
                    << static_cast<WrittenDataType> (tensor[0][2])
                    << static_cast<WrittenDataType> (tensor[1][1])
                    << static_cast<WrittenDataType> (tensor[1][2])
                    << static_cast<WrittenDataType> (tensor[2][2]);
                break;
              }
              case OutputField::Traction:
                xdrWriter << static_cast<WrittenDataType> (dataSource.GetTraction().x)
                    << static_cast<WrittenDataType> (dataSource.GetTraction().y)
                    << static_cast<WrittenDataType> (dataSource.GetTraction().z);
                break;
              case OutputField::TangentialProjectionTraction:
                xdrWriter
                    << static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().x)
                    << static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().y)
                    << static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().z);
                break;
              case OutputField::MpiRank:
                xdrWriter
                    << static_cast<WrittenDataType> (comms.Rank());
                break;
              default:
                // This should never trip. It only occurs when a new OutputField field is added and no
                // implementation is provided for its serialisation.
                assert(false);
            }
          }
        }
      }

      // Actually do the MPI writing.
      outputFile.WriteAt(localDataOffsetIntoFile, buffer);

      // Set the offset to the right place for writing on the next iteration.
      localDataOffsetIntoFile += allCoresWriteLength;
    }

    unsigned LocalPropertyOutput::GetFieldLength(OutputField::FieldType field)
    {
      switch (field)
      {
        case OutputField::Pressure:
        case OutputField::VonMisesStress:
        case OutputField::ShearStress:
        case OutputField::ShearRate:
        case OutputField::MpiRank:
          return 1;
        case OutputField::Velocity:
        case OutputField::Traction:
        case OutputField::TangentialProjectionTraction:
          return 3;
        case OutputField::StressTensor:
          return 6; // We only store the upper triangular part of the symmetric tensor
        default:
          // This should never trip. Only occurs if someone adds a new field and forgets
          // to add to this method.
          assert(false);
          return 0;
      }
    }

    double LocalPropertyOutput::GetOffset(OutputField::FieldType field) const
    {
      switch (field)
      {
        case OutputField::Pressure:
          return REFERENCE_PRESSURE_mmHg;
        default:
          return 0.;
      }
    }
  }
}
