
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_POINTPOINT_IMMEDIATEPOINTPOINT_H
#define HEMELB_NET_MIXINS_POINTPOINT_IMMEDIATEPOINTPOINT_H
#include "net/BaseNet.h"
#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {
    // although ImmediatePointPoint does not use the StoringNet capabilities at all
    // it needs to inherit it
    // so that it becomes the unique final overrider
    class ImmediatePointPoint : public virtual StoringNet
    {

      public:
        ImmediatePointPoint(const MpiCommunicator& comms) :
            BaseNet(comms), StoringNet(comms)
        {
        }
        ~ImmediatePointPoint();

        void WaitPointToPoint();
        // we will *NOT* store the requests, so we must provide RequestSendImpl ourselves.
        virtual void RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type);
        virtual void RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type);
      protected:
        void ReceivePointToPoint(); //PASS
        void SendPointToPoint(); //PASS

    };
  }
}

#endif
