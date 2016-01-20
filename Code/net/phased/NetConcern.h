
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_PHASED_NETCONCERN_H
#define HEMELB_NET_PHASED_NETCONCERN_H

#include "net/phased/Concern.h"
#include "net/phased/steps.h"

namespace hemelb
{
  namespace net
  {
    namespace phased
    {
      class NetConcern : public Concern
      {
        public:
          NetConcern(net::BaseNet & net) :
              net(net)
          {
          }
          bool CallAction(int action)
          {
            switch (static_cast<phased::steps::Step>(action))
            {
              case phased::steps::Send:
                net.Send();
                return true;
              case phased::steps::Receive:
                net.Receive();
                return true;
              case phased::steps::Wait:
                net.Wait();
                return true;

              default:
                return false;
            }
          }
        private:
          net::BaseNet & net;
      };
    }
  }
}

#endif
