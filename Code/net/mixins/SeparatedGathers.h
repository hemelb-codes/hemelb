// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MIXINS_SEPARATEDGATHERS_H
#define HEMELB_NET_MIXINS_SEPARATEDGATHERS_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{

    class SeparatedGathers : public virtual StoringNet
    {
      void ReceiveGathers(){}
      void SendGathers(){}
      void ReceiveGatherVs(){}
      void SendGatherVs(){}
      void WaitGathers();
      void WaitGatherVs();
    };

  }
}
#endif
