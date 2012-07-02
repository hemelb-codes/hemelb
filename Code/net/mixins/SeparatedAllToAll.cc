#include "net/mixins/SeparatedAllToAll.h"
#include <cassert>
namespace hemelb
{
  namespace net
  {
    void SeparatedAllToAll::WaitAllToAll()
    {
      assert(allToAllReceiveProcComms.size()==allToAllSendProcComms.size());
      for (unsigned int i = 0; i < allToAllReceiveProcComms.size(); i++)
      {
        SimpleRequest & sendreq = allToAllSendProcComms[i];
        SimpleRequest & receivereq = allToAllReceiveProcComms[i];
        MPI_Alltoall(sendreq.Pointer,
                     sendreq.Count,
                     sendreq.Type,
                     receivereq.Pointer,
                     receivereq.Count,
                     receivereq.Type,
                     communicator.GetCommunicator());
      }

      allToAllReceiveProcComms.clear();
      allToAllSendProcComms.clear();
    }
  }
}

