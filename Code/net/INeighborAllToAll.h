// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_NET_INEIGHBORALLTOALL_H
#define HEMELB_NET_INEIGHBORALLTOALL_H

#include <vector>
#include <algorithm>
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
    template<class SEND, class RECEIVE = SEND> class INeighborAllToAll
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
        INeighborAllToAll(MpiCommunicator const &comm) :
            comm(comm.Duplicate()), neighbors(comm.GetNeighbors())
        {
        }
        ;
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
        void AddToBuffer(int neighbor, std::vector<Send> const &input);

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

        //! Non-blocking send
        void send();
        //! Blocks until received
        MPI_Status receive();

      protected:
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

      protected:
        //! Returns index of neighbor
        int GetNeighborIndex(int neighbor) const;
    };

    template<class SEND, class RECEIVE> void INeighborAllToAll<SEND, RECEIVE>::send()
    {
      assert(neighbors.size() == 0 or sendBuffer.size() % neighbors.size() == 0);
      assert(comm);
      auto const n = neighbors.size() > 0 ?
        sendBuffer.size() / neighbors.size() :
        0;
      receiveBuffer.resize(sendBuffer.size());

      // Makes sure pointers are valid even if no data
      sendBuffer.reserve(1);
      receiveBuffer.reserve(1);

      // Post message
      HEMELB_MPI_CALL(MPI_Ineighbor_alltoall,
                      ( sendBuffer.data(), n, net::MpiDataType<Send>(), receiveBuffer.data(), n, net::MpiDataType<Receive>(), comm, &request ));
    }

    template<class SEND, class RECEIVE> MPI_Status INeighborAllToAll<SEND, RECEIVE>::receive()
    {
      MPI_Status result;
      HEMELB_MPI_CALL(MPI_Wait, (&request, &result));
      return result;
    }

    template<class SEND, class RECEIVE>
    int INeighborAllToAll<SEND, RECEIVE>::GetNeighborIndex(int neighbor) const
    {
      auto const i_neigh = std::find(neighbors.begin(), neighbors.end(), neighbor);
      assert(i_neigh != neighbors.end());
      return i_neigh - neighbors.begin();
    }

    template<class SEND, class RECEIVE>
    void INeighborAllToAll<SEND, RECEIVE>::AddToBuffer(int neighbor, std::vector<Send> const &input)
    {
      auto const N = input.size();
      sendBuffer.resize(neighbors.size() * N);
      std::copy(input.begin(), input.end(), sendBuffer.begin() + N * GetNeighborIndex(neighbor));
    }

  } /* redblood */
} /* hemelb */
#endif
