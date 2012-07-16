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
namespace hemelb
{
  namespace net
  {
    class Net : public CoalescePointPoint,
                public InterfaceDelegationNet,
                public SeparatedGathers,
                public SeparatedAllToAll
    {
      public:
        Net() :
            BaseNet(), StoringNet()
        {

        }
        Net(topology::Communicator &communicator) :
            BaseNet(communicator), StoringNet()
        {
        }
    };
  }
}
#endif

