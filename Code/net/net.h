#ifndef HEMELB_NET_NET_H
#define HEMELB_NET_NET_H

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

    class Net
    {
      public:
        Net();
        Net(topology::Communicator &communicator);

        ~Net();

        //DTMP: monitoring variables
        long long int BytesSent;
        long long int SyncPointsCounted;

        void Receive();
        void Send();
        void Wait();

        /***
         * Carry out a complete send-receive-wait
         */
        void Dispatch();

        /**
         * Request that iCount entries of type T be included in the send to iToRank,
         * starting at oPointer.
         *
         * @param oPointer Pointer to first element to be sent.
         * @param iCount Number of element to be sent.
         * @param iToRank Rank to send to.
         */
        template<class T>
        void RequestSend(T* oPointer, int iCount, proc_t iToRank);

        /***
         * Request that a signle value be sent
         * @param value Value to send
         * @param toRank Rank to send to
         */
        template<class T>
        void RequestSend(T& value, proc_t toRank);

        /***
         * Request send of the contents of vector specified
         * @param payload Vector whose contents should be sent
         * @param toRank Rank to send to
         */
        template<class T>
        void RequestSend(std::vector<T> &payload, proc_t toRank);

        /**
         * Request that iCount entries of type T be included in the receive from iFromRank
         * starting at oPointer.
         *
         * @param oPointer Pointer to the first element of the array to receive into.
         * @param iCount Number of elements to receive into the array.
         * @param iFromRank Rank to receive from.
         */
        template<class T>
        void RequestReceive(T* oPointer, int iCount, proc_t iFromRank);

        /***
         * Request that a single value be received
         * @param value Value to send
         * @param fromRank Rank to receive from
         */
        template<class T>
        void RequestReceive(T& value, proc_t fromRank);

        /***
         * Request that the vector specified should receive content
         * @param payload Vector to receive content, should already be sized to size of expected message
         * @param fromRank Rank to receive from
         */
        template<class T>
        void RequestReceive(std::vector<T> &payload, proc_t fromRank);

        template<class T>
        void RequestGatherVReceive(T* buffer, int * displacements, int *counts);

        template<class T>
        void RequestGatherVReceive(std::vector<std::vector<T> > &buffer);

        template<class T>
        void RequestGatherReceive(std::vector<T> &buffer);

        template<class T>
        void RequestGatherReceive(T* buffer);

        template<class T>
        void RequestGatherSend(T& value, proc_t toRank);

        template<class T>
        void RequestGatherVSend(std::vector<T> &payload, proc_t toRank);

        template<class T>
        void RequestGatherSend(T* buffer, proc_t toRank);

        template<class T>
        void RequestGatherVSend(T* buffer, int count, proc_t toRank);

      private:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        class BaseRequest
        {
          public:
            void * Pointer;
            int Count;
            MPI_Datatype Type;
            BaseRequest(void *pointer, int count, MPI_Datatype type) :
                Pointer(pointer), Count(count), Type(type)
            {
            }
        };

        class ScalarRequest
        {
          public:
            void * Pointer;
            MPI_Datatype Type;
            ScalarRequest(void *pointer, MPI_Datatype type) :
                Pointer(pointer), Type(type)
            {
            }
        };

        class GatherVReceiveRequest
        {
          public:
            void * Pointer;
            int * Displacements;
            int * Counts;
            MPI_Datatype Type;
            GatherVReceiveRequest(void *pointer,int *displacements, int *counts, MPI_Datatype type) :
            Pointer(pointer),  Displacements(displacements), Counts(counts),Type(type)
            {
            }
        };

        template<class Request>
        class BaseProcComms : public std::vector<Request>
        {
          public:
            MPI_Datatype Type;
        };

        class ProcComms : public BaseProcComms<BaseRequest>
        {
          public:
            void CreateMPIType()
            {
              std::vector<MPI_Aint> displacements(this->size());
              std::vector<int> lengths;
              std::vector<MPI_Datatype> types;

              int lLocation = 0;

              MPI_Aint offset;
              MPI_Get_address(this->front().Pointer, &offset);
              for (typename std::vector<BaseRequest>::iterator it = this->begin(); it != this->end(); ++it)
              {
                MPI_Get_address(it->Pointer, &displacements[lLocation]);
                displacements[lLocation] -= offset;
                ++lLocation;
                lengths.push_back(it->Count);
                types.push_back(it->Type);
              }

              // Create the type and commit it.
              MPI_Type_create_struct(this->size(), &lengths.front(), &displacements.front(), &types.front(), &Type);
              MPI_Type_commit(&Type);
            }
        };

        class GatherProcComms : public BaseProcComms<ScalarRequest>
        {

        };

        class GatherVReceiveProcComms : public BaseProcComms<GatherVReceiveRequest>
        {

        };

        void EnsurePreparedToSendReceive();

        void CreateMPIType(ProcComms *iMetaData);

        void EnsureEnoughRequests(size_t count);

        bool sendReceivePrepped;
        std::map<proc_t, ProcComms> sendProcessorComms;
        std::map<proc_t, ProcComms> receiveProcessorComms;
        std::map<proc_t, ProcComms> gatherVSendProcessorComms;
        GatherVReceiveProcComms gatherVReceiveProcessorComms;
        std::map<proc_t, GatherProcComms> gatherSendProcessorComms;
        GatherProcComms gatherReceiveProcessorComms;

        // Requests and statuses available for general communication within the Net object (both
        // initialisation and during each iteration). Code using these must make sure
        // there are enough available. We do this in a way to minimise the number created
        // on each core, but also to minimise creation / deletion overheads.
        std::vector<MPI_Request> mRequests;
        std::vector<MPI_Status> mStatuses;

        const topology::Communicator &communicator;
    };

  }
}
// include template definitions
#include "net/net.hpp"
#endif // HEMELB_NET_NET_H
