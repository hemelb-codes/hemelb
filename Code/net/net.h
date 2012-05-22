#ifndef HEMELB_NET_NET_H
#define HEMELB_NET_NET_H

#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace net
  {

    class Net
    {
      public:
        Net();
        Net(MPI_Comm commObject);

        ~Net();

        //DTMP: monitoring variables
        long long int BytesSent;
        long long int SyncPointsCounted;

        void Receive();
        void Send();
        void Wait();

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

      private:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        struct ProcComms
        {
          public:
            std::vector<void*> PointerList;
            std::vector<int> LengthList;
            std::vector<MPI_Datatype> TypeList;

            void clear()
            {
              PointerList.clear();
              LengthList.clear();
              TypeList.clear();
            }

            MPI_Datatype Type;
        };

        ProcComms* GetProcComms(proc_t iRank, bool iIsSend);

        template<typename T>
        void AddToList(T* dataToAdd, int dataLength, ProcComms *procCommsObjectToAddTo);

        void EnsurePreparedToSendReceive();

        void CreateMPIType(ProcComms *iMetaData);

        void EnsureEnoughRequests(size_t count);

        bool sendReceivePrepped;

        std::map<proc_t, ProcComms> mSendProcessorComms;
        std::map<proc_t, ProcComms> mReceiveProcessorComms;

        // Requests and statuses available for general communication within the Net object (both
        // initialisation and during each iteration). Code using these must make sure
        // there are enough available. We do this in a way to minimise the number created
        // on each core, but also to minimise creation / deletion overheads.
        std::vector<MPI_Request> mRequests;
        std::vector<MPI_Status> mStatuses;

        MPI_Comm communicator;
    };

  }
}
// include template definitions
#include "net/net.hpp"
#endif // HEMELB_NET_NET_H
