// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

