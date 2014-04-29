// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mixins/alltoall/SeparatedAllToAll.h"
#include <cassert>
namespace hemelb
{
  namespace net
  {
    SeparatedAllToAll::SeparatedAllToAll(const MpiCommunicator& comms) :
        BaseNet(comms), StoringNet(comms)
    {
    }

    void SeparatedAllToAll::WaitAllToAll()
    {
      assert(allToAllReceiveProcComms.size() == allToAllSendProcComms.size());
      for (unsigned int i = 0; i < allToAllReceiveProcComms.size(); i++)
      {
        SimpleRequest & sendreq = allToAllSendProcComms[i];
        SimpleRequest & receivereq = allToAllReceiveProcComms[i];
        HEMELB_MPI_CALL(MPI_Alltoall,
                        (sendreq.Pointer, sendreq.Count, sendreq.Type, receivereq.Pointer, receivereq.Count, receivereq.Type, communicator));
      }

      allToAllReceiveProcComms.clear();
      allToAllSendProcComms.clear();
    }
  }
}

