// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
