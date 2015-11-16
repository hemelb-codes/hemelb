
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_NET_INEIGHBORALLTOALLV_H
#define HEMELB_NET_INEIGHBORALLTOALLV_H

#include <vector>
#include <cassert>
#include <numeric>

#include "net/INeighborAllToAll.h"
#include "net/MpiCommunicator.h"
#include "units.h"

namespace hemelb
{
  namespace net
  {
    //! \brief Simplifies non-blocking neighberhood collectives
    //! \details Owns send and receive buffer. Handles waiting for request to complete.
    //! \code
    //! INeighborAllToAll<element_type> all2all(graph_communicator);
    //! all2all.GetSendBuffer().resize(N * number of neighbors);
    //! all2all.GetSendBuffer()[0] = populate this
    //! all2all.GetSendBuffer()[1] = populate that
    //! ..
    //! all2all.GetSendBuffer()[N * number of neighbors - 1] = populate that
    //! all2all.send() // sends N items to buffer
    //! all2all.receive() // waits until all items are received
    //! all2all.GetReceiveBuffer() // profit!
    //! \endcode
    //! Both SEND and RECEIVE should be types registered with net::MpiDataType;
    template<class SEND, class RECEIVE = SEND>
    class INeighborAllToAllV : protected INeighborAllToAll<SEND, RECEIVE> {
      public:
        //! Sent type
        typedef SEND Send;
        //! Received type
        typedef RECEIVE Receive;
        //! Container holding sent type
        typedef std::vector<Send> SendBuffer;
        //! Container holding received type
        typedef std::vector<Receive> ReceiveBuffer;

        //! Input is a graph communicator
        INeighborAllToAllV(
            MpiCommunicator const &comm,
            std::vector<int> sendCounts,
            std::vector<int> recvCounts)
          : INeighborAllToAll<Send, Receive>(comm.Duplicate()),
            sendCounts(sendCounts), recvCounts(recvCounts)
        {
          sendBuffer.resize(std::accumulate(sendCounts.begin(), sendCounts.end(), 0));
          receiveBuffer.resize(std::accumulate(recvCounts.begin(), recvCounts.end(), 0));
        };
        virtual ~INeighborAllToAllV()
        {
        }

#       define HEMELB_MACRO(NAME, TYPE, CONST)                         \
          TYPE CONST & Get ## NAME()  CONST                            \
          {                                                            \
            return INeighborAllToAll<Send, Receive>::Get ## NAME();    \
          }

          // We don't want exactly the same interface for INeighborAllToAll and INeighborAllToAllV,
          // so redeclare the bits in common and keep the base protected.
          HEMELB_MACRO(SendBuffer, SendBuffer, const);
          HEMELB_MACRO(ReceiveBuffer, ReceiveBuffer, const);
          HEMELB_MACRO(Communicator, MpiCommunicator, const);
#       undef HEMELB_MACRO

#       define HEMELB_MACRO(NAME, VAR)                                \
          void Set ## NAME(std::vector<int> & input)                  \
          {                                                           \
            VAR = input;                                              \
          }                                                           \
          std::vector<int> const & Get ## NAME()  const               \
          {                                                           \
            return VAR;                                               \
          }                                                           \
          int Get ## NAME(proc_t process)                             \
          {                                                           \
            return VAR[GetNeighborIndex(process)];                    \
          }                                                           \
          void Set ## NAME(proc_t process, int count)                 \
          {                                                           \
            VAR[GetNeighborIndex(process)] = count;                   \
          }

          HEMELB_MACRO(SendCounts, sendCounts)
          HEMELB_MACRO(RecvCounts, recvCounts)
#       undef HEMELB_MACRO

        void send();
        MPI_Status receive()
        {
          return INeighborAllToAll<Send, Receive>::receive();
        }

        //! Uses output iterator to fill in send or receive buffer
        template<class ITERATOR> void insertSend(int neighbor, ITERATOR input)
        {
          insert(neighbor, input, sendBuffer, GetSendCounts());
        }
        //! Uses output iterator to fill in send or receive buffer
        void insertSend(int neighbor, std::vector<Send> const & input)
        {
          assert(input.size() == GetSendCounts()[GetNeighborIndex(neighbor)]);
          insert(neighbor, input.begin(), sendBuffer, GetSendCounts());
        }
        //! Uses output iterator to fill in send or receive buffer
        template<class ITERATOR> void insertReceive(int neighbor, ITERATOR input)
        {
          insert(neighbor, input, receiveBuffer, GetRecvCounts());
        }
        //! Uses output iterator to fill in send or receive buffer
        void insertReceive(int neighbor, std::vector<Send> const & input)
        {
          assert(input.size() == GetRecvCounts()[GetNeighborIndex(neighbor)]);
          insert(neighbor, input, receiveBuffer, GetRecvCounts());
        }

      private:
        using INeighborAllToAll<Send, Receive>::sendBuffer;
        using INeighborAllToAll<Send, Receive>::receiveBuffer;
        using INeighborAllToAll<Send, Receive>::request;
        using INeighborAllToAll<Send, Receive>::comm;
        using INeighborAllToAll<Send, Receive>::neighbors;
        using INeighborAllToAll<Send, Receive>::GetNeighborIndex;
        //! Number of objects to send to each proc
        std::vector<int> sendCounts;
        //! Number of objects to receive from each proc
        std::vector<int> recvCounts;

        //! Computes offsets for each neighbor
        static std::vector<int> GetOffsets(std::vector<int> const &counts);

        //! Uses output iterator to fill in send or receive buffer
        template<class ITERATOR, class INPUT>
          void insert(
              int neighbor, ITERATOR input,
              std::vector<INPUT> &container, std::vector<int> const & counts);
    };

    template<class SEND, class RECEIVE>
      std::vector<int>
      INeighborAllToAllV<SEND, RECEIVE>::GetOffsets(std::vector<int> const & counts)
      {
        //! Reserving at least one, so that result.data() is valid.
        std::vector<int> result; result.reserve(1);
        int current = 0;
        for(auto c: counts)
        {
          result.push_back(current);
          current += c;
        }
        return result;
      }

    template<class SEND, class RECEIVE> template<class ITERATOR, class INPUT>
    void INeighborAllToAllV<SEND, RECEIVE>::insert(
        int neighbor, ITERATOR input,
        std::vector<INPUT> &container, std::vector<int> const & counts)
    {
      auto const index = GetNeighborIndex(neighbor);
      assert(counts.size() > index);
      auto const offset = std::accumulate(counts.begin(), counts.begin() + index, 0);
      assert(container.size() >= counts[index] + offset);
      for(int i(offset); i < counts[index] + offset; ++i, ++input)
      {
        container[i] = *input;
      }
    }

    template<class SEND, class RECEIVE>
      void INeighborAllToAllV<SEND, RECEIVE>::send()
      {
        auto const sendOffsets = GetOffsets(sendCounts);
        auto const sendType = net::MpiDataType<Send>();
        auto const recvOffsets = GetOffsets(recvCounts);
        auto const recvType = net::MpiDataType<Receive>();

        // Makes sure buffers are valid even if not communicating anything
        sendBuffer.reserve(1);     sendCounts.reserve(1);
        receiveBuffer.reserve(1);  recvCounts.reserve(1);

        // Post message
        HEMELB_MPI_CALL(
          MPI_Ineighbor_alltoallv,
          (
            sendBuffer.data(), sendCounts.data(), sendOffsets.data(), sendType,
            receiveBuffer.data(), recvCounts.data(), recvOffsets.data(), recvType,
            comm, &request
          )
        );
      }
  } /* redblood */
} /* hemelb */
#endif
