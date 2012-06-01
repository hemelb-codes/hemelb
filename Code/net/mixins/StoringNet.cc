#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {

    void StoringNet::RequestSend(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {

        sendProcessorComms[rank].push_back(BaseRequest(pointer, count, type,rank));
      }
    }

    void StoringNet::RequestReceive(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {
        receiveProcessorComms[rank].push_back(BaseRequest(pointer, count, type,rank));
      }
    }

    /*
     * Blocking gathers are implemented in MPI as a single call for both send/receive
     * But, here we separate send and receive parts, since this interface may one day be used for
     * nonblocking collectives.
     */

    void StoringNet::RequestGatherVSend(void* buffer, int count, proc_t toRank, MPI_Datatype type)
    {
      gatherVSendProcessorComms[toRank].push_back(BaseRequest(buffer, count, type,toRank));
    }

    void StoringNet::RequestGatherReceive(void* buffer, MPI_Datatype type)
    {
      gatherReceiveProcessorComms.push_back(ScalarRequest(buffer, type,0));
    }

    void StoringNet::RequestGatherSend(void* buffer, proc_t toRank, MPI_Datatype type)
    {
      gatherSendProcessorComms[toRank].push_back(ScalarRequest(buffer, type,toRank));
    }

    void StoringNet::RequestGatherVReceive(void* buffer, int * displacements, int *counts, MPI_Datatype type)
    {
      gatherVReceiveProcessorComms.push_back(GatherVReceiveRequest(buffer, displacements, counts, type));
    }
  }
}
