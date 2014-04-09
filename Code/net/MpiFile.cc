//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/MpiFile.h"
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {
    namespace
    {
      void Deleter(MPI_File* fh)
      {
        HEMELB_MPI_CALL(MPI_File_close, (fh));
        delete fh;
      }
    }
    MpiFile::MpiFile() : comm(NULL)
    {

    }

    MpiFile::MpiFile(const MpiCommunicator& parentComm, MPI_File fh) :
        comm(&parentComm)
    {
      filePtr.reset(new MPI_File(fh), Deleter);
    }

    MpiFile MpiFile::Open(const MpiCommunicator& comm, const std::string& filename, int mode,
                          const MPI_Info info)
    {
      MPI_File ans;
      HEMELB_MPI_CALL(
          MPI_File_open,
          (comm, filename.c_str(), mode, info, &ans)
      );
      return MpiFile(comm, ans);
    }


    void MpiFile::Close()
    {
      if (filePtr)
      {
        filePtr.reset();
      }
    }

    MpiFile::operator MPI_File() const
    {
      return *filePtr;
    }

    const MpiCommunicator& MpiFile::GetCommunicator() const
    {
      return *comm;
    }

    void MpiFile::SetView(MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype, const std::string& datarep, MPI_Info info)
    {
      HEMELB_MPI_CALL(
          MPI_File_set_view,
          (*filePtr, disp, etype, filetype, datarep.c_str(), info)
      );
    }
  }
}

