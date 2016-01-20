
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_GATHERS_SEPARATEDGATHERS_H
#define HEMELB_NET_MIXINS_GATHERS_SEPARATEDGATHERS_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{

    class SeparatedGathers : public virtual StoringNet
    {
    public:
      SeparatedGathers(const MpiCommunicator& comms);
    private:
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
