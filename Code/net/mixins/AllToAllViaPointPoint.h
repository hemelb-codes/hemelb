#ifndef HEMELB_NET_MIXINS_ALLTOALLVIAPOINTPOINT_H
#define HEMELB_NET_MIXINS_ALLTOALLVIAPOINTPOINT_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{

    class AllToAllViaPointPoint : public virtual StoringNet
    {
      void ReceiveAllToAll();
      void SendAllToAll();
      void WaitAllToAll(){}
    };
  }
}
#endif
