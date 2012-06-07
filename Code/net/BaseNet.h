#ifndef HEMELB_NET_BASENET_H
#define HEMELB_NET_BASENET_H

#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "mpiInclude.h"
#include "topology/Communicator.h"

namespace hemelb
{
  namespace net
  {

    class BaseNet
    {
      public:
        BaseNet();
        BaseNet(topology::Communicator &communicator);

        virtual ~BaseNet()
        {
        }
        ;

        //DTMP: monitoring variables
        long long int BytesSent;
        long long int SyncPointsCounted;

        void Receive();
        void Send();
        virtual void Wait();

        /***
         * Carry out a complete send-receive-wait
         */
        void Dispatch();

        const topology::Communicator &GetCommunicator() const
        {
          return communicator;
        }
      protected:
        virtual void SendPointToPoint()=0;
        virtual void SendGathers()=0;
        virtual void SendGatherVs()=0;

        virtual void ReceiveGathers()=0;
        virtual void ReceiveGatherVs()=0;
        virtual void ReceivePointToPoint()=0;

        virtual void WaitPointToPoint()=0;
        virtual void WaitGathers()=0;
        virtual void WaitGatherVs()=0;

        virtual void RequestSend(void* pointer, int count, proc_t rank, MPI_Datatype type)=0;
        virtual void RequestReceive(void* pointer, int count, proc_t rank, MPI_Datatype type)=0;


        /*
         * Blocking gathers are implemented in MPI as a single call for both send/receive
         * But, here we separate send and receive parts, since this interface may one day be used for
         * nonblocking collectives.
         */
        virtual void RequestGatherVSend(void* buffer, int count, proc_t toRank, MPI_Datatype type)=0;
        virtual void RequestGatherReceive(void* buffer, MPI_Datatype type)=0;
        virtual void RequestGatherSend(void* buffer, proc_t toRank, MPI_Datatype type)=0;
        virtual void RequestGatherVReceive(void* buffer, int * displacements, int *counts, MPI_Datatype type)=0;

        std::vector<int> & GetDisplacementsBuffer();
        std::vector<int> & GetCountsBuffer();

        const topology::Communicator &communicator;
      private:
        /***
         * Buffers which can be used to store displacements and counts for cleaning up interfaces
         * These will be cleaned up following a Wait/Dispatch
         */
        std::vector<std::vector<int> > displacementsBuffer;
        std::vector<std::vector<int> > countsBuffer;
    };
  }
}
#endif // HEMELB_NET_NET_H
