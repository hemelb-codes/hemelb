//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "geometry/GeometryReaderBase.h"
#include "Exception.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "io/formats/geometry.h"

namespace hemelb
{
  namespace geometry
  {

    GeometryReaderBase::GeometryReaderBase()
    {
    }

    void GeometryReaderBase::OpenFile(const std::string& dataFilePath)
    {
      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      HEMELB_MPI_CALL(MPI_Info_create, (&fileInfo));

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      std::string accessStyle = "access_style";
      std::string accessStyleValue = "sequential";
      std::string buffering = "collective_buffering";
      std::string bufferingValue = "true";

      HEMELB_MPI_CALL(MPI_Info_set, (fileInfo,
              const_cast<char*> (accessStyle.c_str()),
              const_cast<char*> (accessStyleValue.c_str())));
      HEMELB_MPI_CALL(MPI_Info_set, (fileInfo,
              const_cast<char*> (buffering.c_str()),
              const_cast<char*> (bufferingValue.c_str())));

      // Open the file.
      // Stupid C-MPI lack of const-correctness
      try
      {
      HEMELB_MPI_CALL(MPI_File_open,
          (currentComms,
              const_cast<char *>(dataFilePath.c_str()),
              MPI_MODE_RDONLY,
              fileInfo,
              &file)
      );
      }
      catch (Exception& e)
      {
        e << " '" << dataFilePath << "'";
        throw;
      }

      // TODO: Why is there this fflush?
      fflush( NULL);

      // Set the view to the file.
      std::string mode = "native";
      HEMELB_MPI_CALL(
          MPI_File_set_view,
          (file, 0, MPI_CHAR, MPI_CHAR, const_cast<char*> (mode.c_str()), fileInfo)
      );
    }

    /**
     * Read in the section at the beginning of the config file.
     */
    Geometry GeometryReaderBase::ReadPreamble()
    {
      const unsigned preambleBytes = io::formats::geometry::PreambleLength;
      std::vector<char> preambleBuffer = ReadOnAllTasks(preambleBytes);

      // Create an Xdr translator based on the read-in data.
      io::writers::xdr::XdrReader preambleReader = io::writers::xdr::XdrMemReader(&preambleBuffer[0],
                                                                                  preambleBytes);

      unsigned hlbMagicNumber, gmyMagicNumber, version;
      // Read in housekeeping values
      preambleReader.readUnsignedInt(hlbMagicNumber);
      preambleReader.readUnsignedInt(gmyMagicNumber);
      preambleReader.readUnsignedInt(version);

      // Check the value of the HemeLB magic number.
      if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
      {
        throw Exception() << "This file starts with " << hlbMagicNumber
            << "not the HemeLB magic number " << io::formats::HemeLbMagicNumber;
      }

      // Check the value of the geometry file magic number.
      if (gmyMagicNumber != io::formats::geometry::MagicNumber)
      {
        throw Exception() << "This file is not a geometry file: had " << gmyMagicNumber
            << ", not the geometry magic number " << io::formats::geometry::MagicNumber;
      }

      if (version != io::formats::geometry::VersionNumber)
      {
        throw Exception() << "Geometry file version is "<< version
            << ". Currently supported version is " << io::formats::geometry::VersionNumber;
      }

      // Variables we'll read.
      // We use temporary vars here, as they must be the same size as the type in the file
      // regardless of the internal type used.
      unsigned int blocksX, blocksY, blocksZ, blockSize;
      double voxelSize;
      util::Vector3D<double> origin;

      // Read in the values.
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);
      preambleReader.readDouble(voxelSize);
      for (unsigned int i = 0; i < 3; ++i)
      {
        preambleReader.readDouble(origin[i]);
      }

      // Read the padding unsigned int.
      unsigned paddingValue;
      preambleReader.readUnsignedInt(paddingValue);

      return Geometry(util::Vector3D<site_t>(blocksX, blocksY, blocksZ),
                      blockSize,
                      voxelSize,
                      origin);
    }

    std::vector<char> GeometryReaderBase::ReadOnAllTasks(unsigned nBytes)
    {
      std::vector<char> buffer(nBytes);
      if (currentComms.Rank() == HEADER_READING_RANK)
      {
        HEMELB_MPI_CALL(
            MPI_File_read,
            (file, &buffer[0], nBytes, net::MpiDataType(buffer[0]), MPI_STATUS_IGNORE)
        );
      }
      currentComms.Broadcast(buffer, HEADER_READING_RANK);
      return buffer;
    }
  }
}
