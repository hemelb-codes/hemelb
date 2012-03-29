#include <cassert>
#include "extraction/LocalPropertyOutput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace extraction
  {
    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile* outputSpec) :
        dataSource(dataSource), outputSpec(outputSpec)
    {
      // Open the file as write-only, create it if it doesn't exist, don't create if the file
      // already exists.
      MPI_File_open(MPI_COMM_WORLD,
                    const_cast<char*>(outputSpec->filename.c_str()),
                    MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL,
                    MPI_INFO_NULL,
                    &outputFile);

      uint64_t siteCount = 0;

      dataSource.Reset();

      // Count the sites.
      while (dataSource.ReadNext())
      {
        if (outputSpec->geometry->Include(dataSource, dataSource.GetPosition()))
        {
          ++siteCount;
        }
      }

      // Calculate how long local writes need to be.
      //  * 3 longs for the position of a site.
      writeLength = 3 * 8;

      for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
      {
        writeLength += 4 * GetFieldLength(outputSpec->fields[outputNumber].type);
      }

      //  * these are per local site.
      writeLength *= siteCount;

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // The IO proc also writes the iteration number
        writeLength += 8;
      }

      // Calculate the total length written during one iteration.
      MPI_Allreduce(&writeLength, &allCoresWriteLength, 1, MpiDataType<uint64_t>(), MPI_SUM, MPI_COMM_WORLD);

      // Calculate the total number of sites to be written.
      uint64_t allSiteCount = 0;
      MPI_Reduce(&siteCount, &allSiteCount, 1, MpiDataType<uint64_t>(), MPI_SUM, 0, MPI_COMM_WORLD);

      // Write the header information on the IO proc.
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // Three uints for HemeLB magic number, extraction file magic number and
        // version number.
        // Field count (uint)
        // Descriptions (variable strings).
        // Site count (ulong)
        unsigned headerSize = 3 * 4 + 4 + 8;

        for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
        {
          // An extra:
          //  4 bytes are used for the string length in XDR.
          //  4 bytes are used for the number of floats in each field
          size_t stringLength = outputSpec->fields[outputNumber].name.length();
          headerSize += stringLength + 4 + 4;

          // String content is rounded-up to the nearest 4 bytes
          if (stringLength % 4 != 0)
          {
            headerSize += (4 - (stringLength % 4));
          }
        }

        // Create a header buffer, and write it.
        char* headerBuffer = new char[headerSize];
        io::writers::xdr::XdrMemWriter writer(headerBuffer, headerSize);

        writer << (uint32_t) io::formats::HemeLbMagicNumber << (uint32_t) io::formats::Extraction::MagicNumber
            << (uint32_t) io::formats::Extraction::VersionNumber;

        // Write the number of fields and their names and lengths (in float counts)
        writer << (uint32_t) outputSpec->fields.size();
        for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
        {
          writer << outputSpec->fields[outputNumber].name;
        }
        for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
        {
          writer << (uint32_t) GetFieldLength(outputSpec->fields[outputNumber].type);
        }

        // Write the total site count.
        writer << (uint64_t) allSiteCount;

        // Write from the buffer
        MPI_File_write_at(outputFile, 0, headerBuffer, headerSize, MPI_BYTE, MPI_STATUS_IGNORE);

        // And clear it up.
        delete[] headerBuffer;

        // Calculate where each core should start writing - for core 0 this is easy.
        // It passes the value for core 1 to the core.
        localDataOffsetIntoFile = headerSize;

        if (topology::NetworkTopology::Instance()->GetProcessorCount() > 1)
        {
          localDataOffsetIntoFile += writeLength;
          MPI_Send(&localDataOffsetIntoFile, 1, MpiDataType<uint64_t>(), 1, 1, MPI_COMM_WORLD);
          localDataOffsetIntoFile -= writeLength;
        }
      }
      else
      {
        // Receive the writing start position from the previous core.
        MPI_Recv(&localDataOffsetIntoFile,
                 1,
                 MpiDataType<uint64_t>(),
                 topology::NetworkTopology::Instance()->GetLocalRank() - 1,
                 1,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        // Send the next core its start position.
        if (topology::NetworkTopology::Instance()->GetLocalRank()
            != (topology::NetworkTopology::Instance()->GetProcessorCount() - 1))
        {
          localDataOffsetIntoFile += writeLength;
          MPI_Send(&localDataOffsetIntoFile,
                   1,
                   MpiDataType<uint64_t>(),
                   topology::NetworkTopology::Instance()->GetLocalRank() + 1,
                   1,
                   MPI_COMM_WORLD);
          localDataOffsetIntoFile -= writeLength;
        }
      }

      // Create the buffer that we'll write each iteration's data into.
      buffer = new char[writeLength];
    }

    LocalPropertyOutput::~LocalPropertyOutput()
    {
      // Clean up
      MPI_File_close(&outputFile);
      delete[] buffer;
    }

    bool LocalPropertyOutput::ShouldWrite(unsigned long iterationNumber) const
    {
      return ( (iterationNumber % outputSpec->frequency) == 0);
    }

    const PropertyOutputFile* LocalPropertyOutput::GetOutputSpec() const
    {
      return outputSpec;
    }

    void LocalPropertyOutput::Write(unsigned long iterationNumber)
    {
      // Don't write if we shouldn't this iteration.
      if (!ShouldWrite(iterationNumber))
      {
        return;
      }

      // Don't write if this core doesn't do anything.
      if (writeLength <= 0)
      {
        return;
      }

      // Create the buffer.
      io::writers::xdr::XdrMemWriter xdrWriter(buffer, writeLength);

      // Firstly, the IO proc must write the iteration number.
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        xdrWriter << (uint64_t) iterationNumber;
      }

      dataSource.Reset();

      while (dataSource.ReadNext())
      {
        const util::Vector3D<site_t>& position = dataSource.GetPosition();
        if (outputSpec->geometry->Include(dataSource, position))
        {
          // Write the position
          xdrWriter << (uint64_t) position.x << (uint64_t) position.y << (uint64_t) position.z;

          // Write for each field.
          for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
          {
            switch (outputSpec->fields[outputNumber].type)
            {
              case OutputField::Pressure:
                xdrWriter << (float) dataSource.GetPressure();
                break;
              case OutputField::Velocity:
                xdrWriter << (float) dataSource.GetVelocity().x << (float) dataSource.GetVelocity().y
                    << (float) dataSource.GetVelocity().z;
                break;
                // TODO: Work out how to handle the different stresses.
              case OutputField::VonMisesStress:
                xdrWriter << (float) dataSource.GetVonMisesStress();
                break;
              case OutputField::ShearStress:
                xdrWriter << (float) dataSource.GetShearStress();
                break;
              case OutputField::ShearRate:
                xdrWriter << (float) dataSource.GetShearRate();
                break;
            }
          }
        }
      }

      // Actually do the MPI writing.
      MPI_File_write_at(outputFile, localDataOffsetIntoFile, buffer, writeLength, MPI_BYTE, MPI_STATUS_IGNORE);

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
          return 1;
          break;
        case OutputField::Velocity:
          return 3;
          break;
        default:
          // This should never happen. Only occurs if someone adds a new field and forgets
          // to add to this method.
          assert(false);
          return 0;
      }
    }
  }
}
