#ifndef HEMELB_UNITTESTS_NET_NETMOCK_H
#define HEMELB_UNITTESTS_NET_NETMOCK_H
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "mpiInclude.h"
#include "net/net.h"
#include "unittests/net/RecordingNet.h"

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      using namespace hemelb::net;

      class NetMock : public InterfaceDelegationNet,public RecordingNet,
                      public GathersViaPointPoint
      {
        public:
          NetMock(topology::Communicator & communicator) :
              BaseNet(communicator), RecordingNet()
          {
          }

      };
    }
  }
}
#endif // ONCE
