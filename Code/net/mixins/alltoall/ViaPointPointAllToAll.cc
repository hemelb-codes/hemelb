
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/alltoall/ViaPointPointAllToAll.h"
namespace hemelb
{
  namespace net
  {
    ViaPointPointAllToAll::ViaPointPointAllToAll(const MpiCommunicator& comms) :
        BaseNet(comms), StoringNet(comms)
    {
    }
    void ViaPointPointAllToAll::ReceiveAllToAll()
    {

      for (AllToAllProcComms::iterator receivereq = allToAllReceiveProcComms.begin();
          receivereq != allToAllReceiveProcComms.end(); receivereq++)
      {

        int size;
        MPI_Type_size(receivereq->Type, &size);

        for (int source_rank = 0; source_rank < communicator.Size(); source_rank++)
        {

          // The below use of unsigned char is not formally correct (due to the possibility of char not having alignment 1)
          // But we cannot currently see a better solution to avoid compiler warnings from void* arithmetic.
          RequestReceiveImpl(static_cast<unsigned char *>(receivereq->Pointer) + size * source_rank,
                             1,
                             source_rank,
                             receivereq->Type);
        }

      }

      allToAllReceiveProcComms.clear();
    }

    void ViaPointPointAllToAll::SendAllToAll()
    {

      for (AllToAllProcComms::iterator sendreq = allToAllSendProcComms.begin();
          sendreq != allToAllSendProcComms.end(); sendreq++)
      {

        int size;
        MPI_Type_size(sendreq->Type, &size);

        for (int dest_rank = 0; dest_rank < communicator.Size(); dest_rank++)
        {

          RequestSendImpl(static_cast<unsigned char *>(sendreq->Pointer) + size * dest_rank,
                          1,
                          dest_rank,
                          sendreq->Type);
        }

      }

      allToAllSendProcComms.clear();
    }

  }
}
