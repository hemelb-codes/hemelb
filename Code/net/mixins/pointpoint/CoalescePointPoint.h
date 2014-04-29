// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MIXINS_POINTPOINT_COALESCEPOINTPOINT_H
#define HEMELB_NET_MIXINS_POINTPOINT_COALESCEPOINTPOINT_H
#include "net/BaseNet.h"
#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {
    class CoalescePointPoint : public virtual StoringNet
    {

      public:
        CoalescePointPoint(const MpiCommunicator& comms) :
            BaseNet(comms), StoringNet(comms), sendReceivePrepped(false)
        {
        }
        ~CoalescePointPoint();

        void WaitPointToPoint();

      protected:
        void ReceivePointToPoint();
        void SendPointToPoint();

      private:
        void EnsureEnoughRequests(size_t count);
        void EnsurePreparedToSendReceive();
        bool sendReceivePrepped;

        // Requests and statuses available for general communication within the Net object (both
        // initialisation and during each iteration). Code using these must make sure
        // there are enough available. We do this in a way to minimise the number created
        // on each core, but also to minimise creation / deletion overheads.
        std::vector<MPI_Request> requests;
        std::vector<MPI_Status> statuses;
    };
  }
}

#endif
