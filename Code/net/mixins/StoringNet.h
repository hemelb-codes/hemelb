#ifndef HEMELB_NET_MIXINS_STORINGNET_H
#define HEMELB_NET_MIXINS_STORINGNET_H
#include "net/BaseNet.h"
#include "log/Logger.h"
#include "net/ProcComms.h"
namespace hemelb
{
  namespace net
  {
    class StoringNet : public virtual BaseNet
    {
      public:

        void RequestSend(void* pointer, int count, proc_t rank, MPI_Datatype type);
        void RequestReceive(void* pointer, int count, proc_t rank, MPI_Datatype type);

        /*
         * Blocking gathers are implemented in MPI as a single call for both send/receive
         * But, here we separate send and receive parts, since this interface may one day be used for
         * nonblocking collectives.
         */

        void RequestGatherVSend(void* buffer, int count, proc_t toRank, MPI_Datatype type);
        void RequestGatherReceive(void* buffer, MPI_Datatype type);
        void RequestGatherSend(void* buffer, proc_t toRank, MPI_Datatype type);
        void RequestGatherVReceive(void* buffer, int * displacements, int *counts, MPI_Datatype type);
      protected:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */

        std::map<proc_t, ProcComms> sendProcessorComms;
        std::map<proc_t, ProcComms> receiveProcessorComms;
        std::map<proc_t, ProcComms> gatherVSendProcessorComms;
        GatherVReceiveProcComms gatherVReceiveProcessorComms;
        std::map<proc_t, GatherProcComms> gatherSendProcessorComms;
        GatherProcComms gatherReceiveProcessorComms;

    };
  }
}
#endif
