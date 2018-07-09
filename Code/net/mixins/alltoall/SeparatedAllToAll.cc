
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

