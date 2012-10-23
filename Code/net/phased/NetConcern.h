// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
