
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_BASENET_H
#define HEMELB_NET_BASENET_H

#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "comm/Communicator.h"

namespace hemelb
{
  namespace net
  {

    class BaseNet
    {
      public:
        BaseNet(comm::Communicator::ConstPtr communicator);

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

        inline comm::Communicator::ConstPtr GetCommunicator() const
        {
          return communicator;
        }

        inline int Rank() const
        {
          return communicator->Rank();
        }
        inline int Size() const
        {
          return communicator->Size();
        }
      protected:
        virtual void SendPointToPoint()=0;

        virtual void ReceivePointToPoint()=0;

        virtual void WaitPointToPoint()=0;

        // Interfaces exposing MPI_Datatype, not intended for client class use
        virtual void RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)=0;
        virtual void RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)=0;
	
        std::vector<int> & GetDisplacementsBuffer();
        std::vector<int> & GetCountsBuffer();

        comm::Communicator::ConstPtr communicator;
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
