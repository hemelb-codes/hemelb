
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_NET_INEIGHBORALLTOALL_H
#define HEMELB_NET_INEIGHBORALLTOALL_H

#include <vector>
#include <cassert>
#include "net/MpiCommunicator.h"


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
    template<class SEND, class RECEIVE = SEND> class INeighborAllToAll {
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
        INeighborAllToAll(MpiCommunicator const &comm)
          : comm(comm.Duplicate()), neighbors(comm.GetNeighbors())
        {
        };
        virtual ~INeighborAllToAll()
        {
        }

        SendBuffer & GetSendBuffer()
        {
          return sendBuffer;
        }
        SendBuffer const & GetSendBuffer() const
        {
          return sendBuffer;
        }

        //! \brief Adds objects to buffer
        //! \details resizes buffer according to input size (all messages have the same size) and
        //! copies input so it is sent to given neighbor.
        void AddToBuffer(int neighbor, std::vector<Send> const &input)
        {
          auto const i_neigh = std::find(neighbors.begin(), neighbors.end(), neighbor);
          assert(i_neigh != neighbors.end());
          auto const N = input.size();
          auto const i = i_neigh - neighbors.begin();
          sendBuffer.resize(neighbors.size() * N);
          std::copy(input.begin(), input.end(), sendBuffer.begin() + N * i);
        }

        ReceiveBuffer & GetReceiveBuffer()
        {
          return receiveBuffer;
        }
        ReceiveBuffer const & GetReceiveBuffer() const
        {
          return receiveBuffer;
        }

        MpiCommunicator const & GetCommunicator() const
        {
          return comm;
        }

        void send()
        {
          assert(sendBuffer.size() % neighbors.size() == 0);
          auto const n = sendBuffer.size() / neighbors.size();
          receiveBuffer.resize(sendBuffer.size());

          // Post message
          HEMELB_MPI_CALL(
            MPI_Ineighbor_alltoall,
            (
              sendBuffer.data(), n, net::MpiDataType<Send>(),
              receiveBuffer.data(), n, net::MpiDataType<Receive>(),
              comm, &request
            )
          );
        }
        MPI_Status receive()
        {
          MPI_Status result;
          HEMELB_MPI_CALL(MPI_Wait, (&request, &result));
          return result;
        }

      private:
        //! Where information is sent
        SendBuffer sendBuffer;
        //! Where information is receive
        ReceiveBuffer receiveBuffer;
        //! The actual request
        MPI_Request request;
        //! The communicator for this object
        MpiCommunicator comm;
        //! List of neighbors
        std::vector<int> neighbors;
    };
  } /* redblood */
} /* hemelb */
#endif
