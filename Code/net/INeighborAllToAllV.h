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
    class INeighborAllToAllV : protected INeighborAllToAll<SEND, RECEIVE>
    {
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
        INeighborAllToAllV(MpiCommunicator const &comm, std::vector<int> sendCounts,
                           std::vector<int> receiveCounts) :
            INeighborAllToAll<Send, Receive>(comm), sendCounts(sendCounts),
                receiveCounts(receiveCounts)
        {
          sendBuffer.resize(std::accumulate(sendCounts.begin(), sendCounts.end(), 0));
          receiveBuffer.resize(std::accumulate(receiveCounts.begin(), receiveCounts.end(), 0));
        }
        ;
        //! Input is a graph communicator
        INeighborAllToAllV(MpiCommunicator const &comm) :
                INeighborAllToAllV<Send, Receive>(comm,
                                                  std::vector<int>(comm.GetNeighborsCount(), 0),
                                                  std::vector<int>(comm.GetNeighborsCount(), 0))
        {
        }

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
        HEMELB_MACRO(SendBuffer, SendBuffer, const)
        ;HEMELB_MACRO(ReceiveBuffer, ReceiveBuffer, const)
        ;HEMELB_MACRO(Communicator, MpiCommunicator, const)
        ;
#       undef HEMELB_MACRO

#       define HEMELB_MACRO(Name, name)                                                          \
          void Set ## Name ## Counts(std::vector<int> const & input)                             \
          {                                                                                      \
            name ## Counts = input;                                                              \
            auto const N = std::accumulate(name ## Counts.begin(), name ## Counts.end(), 0);     \
            name ## Buffer.resize(N);                                                            \
          }                                                                                      \
          std::vector<int> const & Get ## Name ## Counts()  const                                \
          {                                                                                      \
            return name ## Counts;                                                               \
          }                                                                                      \
          int Get ## Name ## Counts(proc_t process)                                              \
          {                                                                                      \
            assert(int(name ## Counts.size()) > GetNeighborIndex(process));                      \
            return name ## Counts[GetNeighborIndex(process)];                                    \
          }                                                                                      \
          void Set ## Name ## Counts(proc_t process, int count)                                  \
          {                                                                                      \
            assert(int(name ## Counts.size()) > GetNeighborIndex(process));                      \
            name ## Counts[GetNeighborIndex(process)] = count;                                   \
            auto const N = std::accumulate(name ## Counts.begin(), name ## Counts.end(), 0);     \
            name ## Buffer.resize(N);                                                            \
          }                                                                                      \
          /** Gets nth object from/for a given neighbor  **/                                     \
          Name const & Get ## Name(int neighbor, int i) const                                    \
          {                                                                                      \
            auto const index = GetNeighborIndex(neighbor);                                       \
            assert(int(name ## Counts.size()) > index);                                          \
            auto const offset = std::accumulate(                                                 \
                name ## Counts.begin(), name ## Counts.begin() + index, 0) + i;                  \
            assert(int(name ## Buffer.size()) > offset);                                         \
            return name ## Buffer[offset];                                                       \
          }                                                                                      \
          /** Sets specific send object **/                                                      \
          void Set ## Name(int neighbor, Name const &input, int i)                               \
          {                                                                                      \
            SetInternal(neighbor, input, name ## Buffer, name ## Counts, i);                     \
          }                                                                                      \
          /** Uses output iterator to fill in send or receive buffer **/                         \
          template<class ITERATOR> void insert ## Name(int neighbor, ITERATOR input)             \
          {                                                                                      \
            insert(neighbor, input, name ## Buffer, name ## Counts);                             \
          }                                                                                      \
          /** Uses output iterator to fill in send or receive buffer **/                         \
          void insert ## Name(int neighbor, std::vector<Name> const & input)                     \
          {                                                                                      \
            assert(int(input.size()) == name ## Counts[GetNeighborIndex(neighbor)]);             \
            insert(neighbor, input.begin(), name ## Buffer, name ## Counts);                     \
          }                                                                                      \
          /** Fill buffer with specific input **/                                                \
          template<class INPUT> void fill ## Name(INPUT input)                                   \
          {                                                                                      \
            std::fill(name ## Buffer.begin(), name ## Buffer.end(), input);                      \
          }

        HEMELB_MACRO(Send, send)
        ;HEMELB_MACRO(Receive, receive)
        ;
#       undef HEMELB_MACRO

        void send();
        MPI_Status receive()
        {
          return INeighborAllToAll<Send, Receive>::receive();
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
        std::vector<int> receiveCounts;

        //! Computes offsets for each neighbor
        static std::vector<int> GetOffsets(std::vector<int> const &counts);

        //! Uses output iterator to fill in send or receive buffer
        template<class ITERATOR, class INPUT>
        void insert(int neighbor, ITERATOR input, std::vector<INPUT> &container,
                    std::vector<int> const & counts);

        //! Sets give object
        template<class INPUT>
        void SetInternal(int neighbor, INPUT const &input, std::vector<INPUT> &container,
                         std::vector<int> const & counts, int i);
    };

    template<class SEND, class RECEIVE>
    std::vector<int> INeighborAllToAllV<SEND, RECEIVE>::GetOffsets(std::vector<int> const & counts)
    {
      //! Reserving at least one, so that result.data() is valid.
      std::vector<int> result;
      result.reserve(1);
      int current = 0;
      for (auto c : counts)
      {
        result.push_back(current);
        current += c;
      }
      return result;
    }

    template<class SEND, class RECEIVE> template<class ITERATOR, class INPUT>
    void INeighborAllToAllV<SEND, RECEIVE>::insert(int neighbor, ITERATOR input,
                                                   std::vector<INPUT> &container,
                                                   std::vector<int> const & counts)
    {
      auto const index = GetNeighborIndex(neighbor);
      assert(int(counts.size()) > index);
      auto const offset = std::accumulate(counts.begin(), counts.begin() + index, 0);
      assert(int(container.size()) >= counts[index] + offset);
      for (int i(offset); i < counts[index] + offset; ++i, ++input)
      {
        container[i] = *input;
      }
    }

    template<class SEND, class RECEIVE> template<class INPUT>
    void INeighborAllToAllV<SEND, RECEIVE>::SetInternal(int neighbor, INPUT const &input,
                                                        std::vector<INPUT> &container,
                                                        std::vector<int> const & counts, int i)
    {
      auto const index = GetNeighborIndex(neighbor);
      assert(int(counts.size()) > index);
      auto const offset = std::accumulate(counts.begin(), counts.begin() + index, 0) + i;
      assert(int(container.size()) > offset);
      container[offset] = input;
    }

    template<class SEND, class RECEIVE>
    void INeighborAllToAllV<SEND, RECEIVE>::send()
    {
      auto const sendOffsets = GetOffsets(sendCounts);
      auto const sendType = net::MpiDataType<Send>();
      auto const recvOffsets = GetOffsets(receiveCounts);
      auto const recvType = net::MpiDataType<Receive>();

      // Makes sure buffers are valid even if not communicating anything
      sendBuffer.reserve(1);
      sendCounts.reserve(1);
      receiveBuffer.reserve(1);
      receiveCounts.reserve(1);

      // Post message
      HEMELB_MPI_CALL(MPI_Ineighbor_alltoallv,
                      ( sendBuffer.data(), sendCounts.data(), sendOffsets.data(), sendType, receiveBuffer.data(), receiveCounts.data(), recvOffsets.data(), recvType, comm, &request ));
    }
  } /* redblood */
} /* hemelb */
#endif
