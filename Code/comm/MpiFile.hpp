
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_MPIFILE_HPP
#define HEMELB_COMM_MPIFILE_HPP

#include "comm/MpiFile.h"

namespace hemelb
{
  namespace comm
  {

    template<typename T>
    void MpiFile::Read(std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_read,
          (*filePtr, buffer.data(), buffer.size(), MpiDataType<T>(), stat)
      );
    }
    template<typename T>
    void MpiFile::ReadAt(MPI_Offset offset, std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_read_at,
          (*filePtr, offset, buffer.data(), buffer.size(), MpiDataType<T>(), stat)
      );
    }

    template<typename T>
    void MpiFile::Write(const std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write,
          (*filePtr, buffer.data(), buffer.size(), MpiDataType<T>(), stat)
      );
    }
    template<typename T>
    void MpiFile::WriteAt(MPI_Offset offset, const std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write_at,
          (*filePtr, offset, buffer.data(), buffer.size(), MpiDataType<T>(), stat)
      );

    }

  }
}

#endif
