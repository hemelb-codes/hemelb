
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <mpi.h>

#include "comm/MpiEnvironment.h"
#include "comm/MpiError.h"
#include "comm/MpiCommunicator.h"

namespace hemelb
{
  namespace comm
  {

    MpiEnvironment::MpiEnvironment(int& argc, char**& argv) :
      doesOwnMpi(false)
    {
      if (!Initialized())
      {
        HEMELB_MPI_CALL(MPI_Init, (&argc, &argv));
        HEMELB_MPI_CALL(MPI_Comm_set_errhandler, (MPI_COMM_WORLD, MPI_ERRORS_RETURN));
        doesOwnMpi = true;
      }
    }

    MpiEnvironment::~MpiEnvironment()
    {
      if (doesOwnMpi)
      {
        HEMELB_MPI_CALL_NOTHROW(MPI_Finalize, ());
      }
    }

    Communicator::Ptr MpiEnvironment::World()
    {
      return std::make_shared<MpiCommunicator>(MPI_COMM_WORLD, false);
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
      World()->Abort(errorCode);
    }

  }
}
