//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_GEOMETRY_GEOMETRYREADERBASE_H
#define HEMELB_GEOMETRY_GEOMETRYREADERBASE_H

#include <string>
#include <vector>

#include "net/mpi.h"
#include "geometry/Geometry.h"

namespace hemelb
{
  namespace geometry
  {
    class GeometryReaderBase
    {
      public:
        GeometryReaderBase();
      protected:
        /**
         * Open the GMY file ready for reading. Note must set currentComms before calling this.
         * @param dataFilePath
         */
        void OpenFile(const std::string& dataFilePath);

        Geometry ReadPreamble();
        /**
         * Read from the file into a buffer. We read this on a single core then broadcast it.
         * This has proven to be more efficient than reading in on every core (even using a collective
         * read).
         *
         * Note this allocates memory and returns the pointer to you.
         *
         * @param nBytes
         * @return
         */
        std::vector<char> ReadOnAllTasks(unsigned nBytes);

        //! The rank which reads in the header information.
        static const proc_t HEADER_READING_RANK = 0;

        //! File accessed to read in the geometry data.
        MPI_File file;
        //! Information about the file, to give cues and hints to MPI.
        MPI_Info fileInfo;
        // TODO: This was never a good plan, better code design will avoid the need for it.
        net::MpiCommunicator currentComms; //! The communicator currently in use.
    };

  }
}
#endif // HEMELB_GEOMETRY_GEOMETRYREADERBASE_H
