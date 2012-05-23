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
        const topology::Communicator &communicator;

        void Receive();
        void Send();
        virtual void Wait();

        /***
         * Carry out a complete send-receive-wait
         */
        void Dispatch();

        /*****************************************************************************
         * The interface for net also includes a large number of templated methods which should be defined in the implementing classes
         * We cannot formally write these here as virtual templates, because of a slight suckyness in C++
         * But I list them here for the record.

         template<class T> void RequestSend(T* oPointer, int iCount, proc_t iToRank);
         template<class T> void RequestSend(T& value, proc_t toRank);
         template<class T> void RequestSend(std::vector<T> &payload, proc_t toRank);
         template<class T> void RequestReceive(T* oPointer, int iCount, proc_t iFromRank);
         template<class T> void RequestReceive(T& value, proc_t fromRank);
         template<class T> void RequestReceive(std::vector<T> &payload, proc_t fromRank);

         template<class T> void RequestGatherVReceive(T* buffer, int * displacements, int *counts);
         template<class T> void RequestGatherVReceive(std::vector<std::vector<T> > &buffer);
         template<class T> void RequestGatherVSend(T* buffer, int count, proc_t toRank);
         template<class T> void RequestGatherVSend(T* buffer, int count, proc_t toRank);

         template<class T> void RequestGatherReceive(std::vector<T> &buffer);
         template<class T> void RequestGatherReceive(T* buffer);
         template<class T> void RequestGatherSend(T& value, proc_t toRank);
         template<class T> void RequestGatherSend(T* buffer, proc_t toRank);

         */

      protected:
        virtual void SendPointToPoint()=0;
        virtual void SendGathers()=0;
        virtual void SendGatherVs()=0;
        virtual void ReceiveGathers()=0;
        virtual void ReceiveGatherVs()=0;
        virtual void ReceivePointToPoint()=0;
        virtual void EnsurePreparedToSendReceive()=0;
    };
  }
}
#endif // HEMELB_NET_NET_H
