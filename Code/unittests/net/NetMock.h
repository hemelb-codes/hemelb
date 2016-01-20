
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_NET_NETMOCK_H
#define HEMELB_UNITTESTS_NET_NETMOCK_H
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "net/mpi.h"
#include "net/net.h"
#include "unittests/net/RecordingNet.h"

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      using namespace hemelb::net;

      class NetMock : public InterfaceDelegationNet,
                      public RecordingNet,
                      public ViaPointPointGathers,
                      public ViaPointPointAllToAll
      {
        public:
          NetMock(net::MpiCommunicator & communicator) :
              BaseNet(communicator), StoringNet(communicator), InterfaceDelegationNet(communicator),
                  RecordingNet(communicator), ViaPointPointGathers(communicator),
                  ViaPointPointAllToAll(communicator)
          {
          }

      };
    }
  }
}
#endif // ONCE
