
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_POINTPOINT_SEPARATEDPOINTPOINT_H
#define HEMELB_NET_MIXINS_POINTPOINT_SEPARATEDPOINTPOINT_H
#include "net/BaseNet.h"
#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {
    class SeparatedPointPoint : public virtual StoringNet
    {

      public:
        SeparatedPointPoint(const MpiCommunicator& comms) :
            BaseNet(comms), StoringNet(comms), sendReceivePrepped(false), count_sends(0),
                count_receives(0)
        {
        }
        ~SeparatedPointPoint();

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
        unsigned int count_sends;
        unsigned int count_receives;
    };
  }
}

#endif
