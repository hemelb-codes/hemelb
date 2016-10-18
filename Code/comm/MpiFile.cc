
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/MpiFile.h"
#include "net/MpiCommunicator.h"
#include "net/MpiConstness.h"

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
          (comm, MpiConstCast(filename.c_str()), mode, info, &ans)
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

    MPI_Offset MpiFile::GetSize() const {
      MPI_Offset ans;
      HEMELB_MPI_CALL(
		      MPI_File_get_size,
		      (*filePtr, &ans)
		      );
      return ans;
    }

    const MpiCommunicator& MpiFile::GetCommunicator() const
    {
      return *comm;
    }

    void MpiFile::SetView(MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype, const std::string& datarep, MPI_Info info)
    {
      HEMELB_MPI_CALL(
          MPI_File_set_view,
          (*filePtr, disp, etype, filetype, MpiConstCast(datarep.c_str()), info)
      );
    }
  }
}

