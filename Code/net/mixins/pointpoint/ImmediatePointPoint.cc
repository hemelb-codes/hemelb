
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/pointpoint/ImmediatePointPoint.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace net
  {
    void ImmediatePointPoint::RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      HEMELB_MPI_CALL(
          MPI_Ssend,
          (pointer, count, type, rank, 10, communicator)
      );
    }
    void ImmediatePointPoint::RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      HEMELB_MPI_CALL(
          MPI_Recv,
          (pointer, count, type, rank, 10, communicator, MPI_STATUS_IGNORE)
      );
    }

    void ImmediatePointPoint::ReceivePointToPoint()
    {

    }

    void ImmediatePointPoint::SendPointToPoint()
    {

    }

    /*!
     Free the allocated data.
     */
    ImmediatePointPoint::~ImmediatePointPoint()
    {
    }

    void ImmediatePointPoint::WaitPointToPoint()
    {

    }
  }
}
