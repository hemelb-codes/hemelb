
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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
    void MpiFile::Write(const std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write,
          (*filePtr, MpiConstCast(&buffer[0]), buffer.size(), MpiDataType<T>(), stat)
      );
    }
    template<typename T>
    void MpiFile::WriteAt(MPI_Offset offset, const std::vector<T>& buffer, MPI_Status* stat)
    {
      HEMELB_MPI_CALL(
          MPI_File_write_at,
          (*filePtr, offset, MpiConstCast(&buffer[0]), buffer.size(), MpiDataType<T>(), stat)
      );

    }

  }
}

#endif
