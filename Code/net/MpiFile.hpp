//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPIFILE_HPP
#define HEMELB_NET_MPIFILE_HPP

#include "net/MpiFile.h"

namespace hemelb
{
  namespace net
  {

    template<typename T>
    void MpiFile::Read(std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_read,
          (*filePtr, &buffer[0], buffer.size(), MpiDataType<T>(), stat)
      );
    }
    template<typename T>
    void MpiFile::ReadAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_read_at,
          (*filePtr, offset, &buffer[0], buffer.size(), MpiDataType<T>(), stat)
      );
    }

    template<typename T>
    void MpiFile::Write(std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write,
          (*filePtr, &buffer[0], buffer.size(), MpiDataType<T>(), stat)
      );
    }
    template<typename T>
    void MpiFile::WriteAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write_at,
          (*filePtr, offset, &buffer[0], buffer.size(), MpiDataType<T>(), stat)
      );

    }

  }
}

#endif
