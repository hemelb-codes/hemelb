// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mixins/AllToAllViaPointPoint.h"
namespace hemelb
{
  namespace net
  {
    void AllToAllViaPointPoint::ReceiveAllToAll()
    {

      for (AllToAllProcComms::iterator receivereq = allToAllReceiveProcComms.begin();
          receivereq != allToAllReceiveProcComms.end(); receivereq++)
      {

        int size;
        MPI_Type_size(receivereq->Type, &size);

        for (int source_rank = 0; source_rank < communicator.GetSize(); source_rank++)
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

    void AllToAllViaPointPoint::SendAllToAll()
    {

      for (AllToAllProcComms::iterator sendreq = allToAllSendProcComms.begin(); sendreq != allToAllSendProcComms.end();
          sendreq++)
      {

        int size;
        MPI_Type_size(sendreq->Type, &size);

        for (int dest_rank = 0; dest_rank < communicator.GetSize(); dest_rank++)
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
