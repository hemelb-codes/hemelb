
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_NET_H
#define HEMELB_NET_NET_H

#include "net/BaseNet.h"
#include "net/mixins/mixins.h"
#include "net/BuildInfo.h"
namespace hemelb
{
  namespace net
  {
    class Net : public PointPointImpl,
                public InterfaceDelegationNet,
                public AllToAllImpl,
                public GathersImpl
    {
      public:
        Net(const MpiCommunicator &communicator) :
            BaseNet(communicator), StoringNet(communicator), PointPointImpl(communicator),
                InterfaceDelegationNet(communicator), AllToAllImpl(communicator),
                GathersImpl(communicator)
        {
        }
    };
  }
}
#endif

