//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#include <mpi.h>

#include "net/MpiEnvironment.h"
#include "net/MpiError.h"

namespace hemelb
{
  namespace net
  {

    MpiEnvironment::MpiEnvironment(int& argc, char**& argv) :
      doesOwnMpi(false)
    {
      if (!Initialized())
      {
        HEMELB_MPI_CALL(MPI_Init, (&argc, &argv));
        HEMELB_MPI_CALL(MPI_Errhandler_set, (MPI_COMM_WORLD, MPI_ERRORS_RETURN));
        doesOwnMpi = true;
      }
    }

    MpiEnvironment::~MpiEnvironment()
    {
      if (doesOwnMpi)
      {
        HEMELB_MPI_CALL(MPI_Finalize, ());
      }
    }

    bool MpiEnvironment::Initialized()
    {
      int flag;
      MPI_Initialized(&flag);
      return flag
        ? true
        : false;
    }

    bool MpiEnvironment::Finalized()
    {
      int flag;
      MPI_Finalized(&flag);
      return flag
        ? true
        : false;
    }

    void MpiEnvironment::Abort(int errorCode)
    {
      MPI_Abort(MPI_COMM_WORLD, errorCode);
    }

  }
}
