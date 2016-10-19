
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {
    StoringNet::StoringNet(comm::Communicator::ConstPtr comms) :
        BaseNet(comms)
    {
    }
    void StoringNet::RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {
        sendProcessorComms[rank].push_back(SimpleRequest(pointer, count, type, rank));
      }
    }

    void StoringNet::RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {
        receiveProcessorComms[rank].push_back(SimpleRequest(pointer, count, type, rank));
      }
    }

  }
}
