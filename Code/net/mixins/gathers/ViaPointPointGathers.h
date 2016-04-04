
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MIXINS_GATHERS_VIAPOINTPOINTGATHERS_H
#define HEMELB_NET_MIXINS_GATHERS_VIAPOINTPOINTGATHERS_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{
    /***
     * Reimplement gathers via point-point calls
     * This code is not robust for working with gathers of complex defined datatypes.
     * It is intended for testing and performance measurement use, and should not be used in production.
     */
    class ViaPointPointGathers : public virtual StoringNet
    {
    public:
      ViaPointPointGathers(const MpiCommunicator& comms);
    private:
      void ReceiveGathers();
      void SendGathers();
      void ReceiveGatherVs();
      void SendGatherVs();
      void WaitGathers(){}
      void WaitGatherVs(){}
    };
  }
}
#endif
