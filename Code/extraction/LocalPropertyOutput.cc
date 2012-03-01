#include "extraction/LocalPropertyOutput.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace extraction
  {
    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile& outputSpec) :
      dataSource(dataSource), outputSpec(outputSpec)
    {
      // Open the file as write-only, create it if it doesn't exist, don't create if the file
      // already exists.
      MPI_File_open(MPI_COMM_WORLD, const_cast<char*> (outputSpec.filename.c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE
          | MPI_MODE_EXCL, MPI_INFO_NULL, &outputFile);

      buffer = new char[writeLength];

      uint64_t siteCount = 0;

      dataSource.Reset();

      util::Vector3D<site_t> position;
      util::Vector3D<float> velocity;
      float pressure, stress;

      while (dataSource.ReadNext(position, pressure, velocity, stress))
      {
        if (outputSpec.geometry->Include(dataSource, position))
        {
          ++siteCount;
        }
      }

      writeLength = 0;

      for (unsigned outputNumber = 0; outputNumber < outputSpec.fields.size(); ++outputNumber)
      {
        switch (outputSpec.fields[outputNumber].type)
        {
          case OutputField::Pressure:
            writeLength += 4;
            break;
          case OutputField::Velocity:
            writeLength += 12;
            break;
            // TODO: Work out how to handle the different stresses.
          case OutputField::VonMisesStress:
          case OutputField::ShearStress:
            writeLength += 4;
            break;
        }
      }

      writeLength *= siteCount;

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // Iteration number
        writeLength += 8;
      }

      MPI_Allreduce(&writeLength, &allCoresWriteLength, 1, MpiDataType<unsigned> (), MPI_SUM, MPI_COMM_WORLD);

      uint64_t allSiteCount = 0;
      MPI_Reduce(&siteCount, &allSiteCount, 1, MpiDataType<uint64_t> (), MPI_SUM, 0, MPI_COMM_WORLD);

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // Field count (uint)
        // Descriptions (variable strings).
        // Site count (ulong)
        unsigned headerSize = 4 + 8;

        for (unsigned outputNumber = 0; outputNumber < outputSpec.fields.size(); ++outputNumber)
        {
          headerSize += outputSpec.fields[outputNumber].name.length() + 1;
        }

        char* headerBuffer = new char[headerSize];
        io::writers::xdr::XdrMemWriter writer(headerBuffer, headerSize);

        writer << outputSpec.fields.size();
        for (unsigned outputNumber = 0; outputNumber < outputSpec.fields.size(); ++outputNumber)
        {
          writer << outputSpec.fields[outputNumber].name;
        }
        writer << (uint64_t) allSiteCount;

        MPI_File_write_at(outputFile, 0, headerBuffer, headerSize, MPI_BYTE, MPI_STATUS_IGNORE);

        delete[] headerBuffer;

        localDataOffsetIntoFile = headerSize;

        localDataOffsetIntoFile += writeLength;
        MPI_Send(&localDataOffsetIntoFile, 1, MpiDataType<MPI_Offset> (), 1, 1, MPI_COMM_WORLD);
        localDataOffsetIntoFile -= writeLength;
      }
      else
      {
        MPI_Recv(&localDataOffsetIntoFile,
                 1,
                 MpiDataType<MPI_Offset> (),
                 topology::NetworkTopology::Instance()->GetLocalRank() - 1,
                 1,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        if (topology::NetworkTopology::Instance()->GetLocalRank()
            != topology::NetworkTopology::Instance()->GetProcessorCount() - 1)
        {
          localDataOffsetIntoFile += writeLength;
          MPI_Send(&localDataOffsetIntoFile,
                   1,
                   MpiDataType<MPI_Offset> (),
                   topology::NetworkTopology::Instance()->GetLocalRank() + 1,
                   1,
                   MPI_COMM_WORLD);
          localDataOffsetIntoFile -= writeLength;
        }
      }
    }

    LocalPropertyOutput::~LocalPropertyOutput()
    {
      MPI_File_close(&outputFile);
      delete[] buffer;
    }

    void LocalPropertyOutput::Write(unsigned long iterationNumber)
    {
      if ( (iterationNumber % outputSpec.frequency) != 0)
      {
        return;
      }

      if (writeLength <= 0)
      {
        return;
      }

      io::writers::xdr::XdrMemWriter xdrWriter(buffer, writeLength);

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        xdrWriter << (uint64_t) iterationNumber;
      }

      util::Vector3D<site_t> position;
      util::Vector3D<float> velocity;
      float pressure, stress;

      dataSource.Reset();

      while (dataSource.ReadNext(position, pressure, velocity, stress))
      {
        if (outputSpec.geometry->Include(dataSource, position))
        {
          xdrWriter << position.x << position.y << position.z;

          for (unsigned outputNumber = 0; outputNumber < outputSpec.fields.size(); ++outputNumber)
          {
            switch (outputSpec.fields[outputNumber].type)
            {
              case OutputField::Pressure:
                xdrWriter << pressure;
                break;
              case OutputField::Velocity:
                xdrWriter << velocity.x << velocity.y << velocity.z;
                break;
                // TODO: Work out how to handle the different stresses.
              case OutputField::VonMisesStress:
              case OutputField::ShearStress:
                xdrWriter << stress;
                break;
            }
          }
        }
      }

      MPI_File_write_at(outputFile, localDataOffsetIntoFile, buffer, writeLength, MPI_BYTE, MPI_STATUS_IGNORE);

      localDataOffsetIntoFile += allCoresWriteLength;
    }
  }
}
