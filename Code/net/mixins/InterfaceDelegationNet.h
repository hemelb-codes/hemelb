// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_H
#define HEMELB_NET_MIXINS_INTERFACEDELEGATIONNET_H
namespace hemelb
{
  namespace net
  {
    /***
     * Define the external template-based interface to be used by client classes.
     * We define a series of templates with nice C++ style interfaces which delegate to C-style template interfaces
     * And then, delegate the templated C-style interfaces to nontemplated interfaces taking an MPI datatype.
     */
    class InterfaceDelegationNet : public virtual BaseNet
    {
      public:
        InterfaceDelegationNet(const MpiCommunicator& comms) :
            BaseNet(comms)
        {
        }

        template<class T>
        void RequestSendV(std::vector<T> &payload, proc_t toRank)
        {
          RequestSend(&payload[0], payload.size(), toRank);
        }

        template<class T>
        void RequestSendR(T& value, proc_t toRank)
        {
          RequestSend(&value, 1, toRank);
        }

        template<class T>
        void RequestReceiveR(T& value, proc_t fromRank)
        {
          RequestReceive(&value, 1, fromRank);
        }

        template<class T>
        void RequestReceiveV(std::vector<T> &payload, proc_t toRank)
        {
          RequestReceive(&payload[0], payload.size(), toRank);
        }

        template<class T>
        void RequestGatherVReceive(std::vector<T>& bigBuffer, const std::vector<int>& countsIn)
        {
          // It may seem inefficient to go through the list twice (it is), but it avoids any nastiness
          // in 32 / 64 bit problems, because we can take our displacement as simply the difference
          // between two addresses (by allocating the bigBuffer before the address calculation).

          int totalCount = 0;
          for (std::vector<int>::const_iterator count_iterator = countsIn.begin();
              count_iterator != countsIn.end(); count_iterator++)
          {
            totalCount += *count_iterator;
          }

          // Allocate the large buffer we'll receive all data into.
          bigBuffer.resize(totalCount);
          bigBuffer.reserve(1);

          std::vector<int> & displacements = this->GetDisplacementsBuffer();
          std::vector<int> &counts = this->GetCountsBuffer();

          // Now store the displacement and count for each sending proc.
          int countSoFar = 0;

          for (std::vector<int>::const_iterator count_iterator = countsIn.begin();
              count_iterator != countsIn.end(); count_iterator++)
          {
            int nextCount = *count_iterator;

            counts.push_back(nextCount);
            displacements.push_back(&bigBuffer[countSoFar] - &bigBuffer[0]);

            countSoFar += nextCount;
          }

          // And store the pointer and culminating arrays.
          RequestGatherVReceive(&bigBuffer.front(), &displacements.front(), &counts.front());
        }

        template<class T>
        void RequestGatherReceive(std::vector<T> &buffer)
        {
          // Ensure vector has some underlying array, even if it's unused.
          buffer.reserve(1);
          RequestGatherReceive(&buffer.front());
        }

        template<class T>
        void RequestGatherSend(T& value, proc_t toRank)
        {
          RequestGatherSend(&value, toRank);
        }

        template<class T>
        void RequestGatherVSend(std::vector<T> &payload, proc_t toRank)
        {
          RequestGatherVSend(&payload.front(), payload.size(), toRank);
        }

        /***
         * This is for a scalar all to all
         * @param buffer vector with length same as communicator size
         */
        template<class T>
        void RequestAllToAllSend(std::vector<T> &buffer)
        {
          RequestAllToAllSend(&buffer.front(), 1);
        }
        /***
         * This is for a scalar all to all
         * @param buffer vector with length same as communicator size
         */
        template<class T>
        void RequestAllToAllReceive(std::vector<T> &buffer)
        {
          RequestAllToAllReceive(&buffer.front(), 1);
        }

        template<class T>
        void RequestSend(T* pointer, int count, proc_t rank)
        {
          RequestSendImpl(pointer, count, rank, MpiDataType<T>());
        }

        template<class T>
        void RequestReceive(T* pointer, int count, proc_t rank)
        {
          RequestReceiveImpl(pointer, count, rank, MpiDataType<T>());
        }

        /*
         * Blocking gathers are implemented in MPI as a single call for both send/receive
         * But, here we separate send and receive parts, since this interface may one day be used for
         * nonblocking collectives.
         */

        template<class T>
        void RequestGatherVSend(T* buffer, int count, proc_t toRank)
        {
          RequestGatherVSendImpl(buffer, count, toRank, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherReceive(T* buffer)
        {
          RequestGatherReceiveImpl(buffer, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherSend(T* buffer, proc_t toRank)
        {
          RequestGatherSendImpl(buffer, toRank, MpiDataType<T>());
        }

        template<class T>
        void RequestGatherVReceive(T* buffer, int * displacements, int *counts)
        {
          RequestGatherVReceiveImpl(buffer, displacements, counts, MpiDataType<T>());
        }

        template<class T>
        void RequestAllToAllSend(T* buffer, int count)
        {
          RequestAllToAllSendImpl(buffer, count, MpiDataType<T>());
        }
        template<class T>
        void RequestAllToAllReceive(T* buffer, int count)
        {
          RequestAllToAllReceiveImpl(buffer, count, MpiDataType<T>());
        }

    }
    ;
  }
}
#endif
