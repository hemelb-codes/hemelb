
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_ALLTOALL_VIAPOINTPOINTALLTOALL_H
#define HEMELB_NET_MIXINS_ALLTOALL_VIAPOINTPOINTALLTOALL_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{

    class ViaPointPointAllToAll : public virtual StoringNet
    {
    public:
      ViaPointPointAllToAll(const MpiCommunicator& comms);
    private:
      void ReceiveAllToAll();
      void SendAllToAll();
      void WaitAllToAll(){}
    };
  }
}
#endif
