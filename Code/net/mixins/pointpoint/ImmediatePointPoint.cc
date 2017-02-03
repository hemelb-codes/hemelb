// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mixins/pointpoint/ImmediatePointPoint.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace net
  {
    void ImmediatePointPoint::RequestSendImpl(void* pointer, int count, proc_t rank,
                                              MPI_Datatype type)
    {
      HEMELB_MPI_CALL(MPI_Ssend, (pointer, count, type, rank, 10, communicator));
    }
    void ImmediatePointPoint::RequestReceiveImpl(void* pointer, int count, proc_t rank,
                                                 MPI_Datatype type)
    {
      HEMELB_MPI_CALL(MPI_Recv, (pointer, count, type, rank, 10, communicator, MPI_STATUS_IGNORE));
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
