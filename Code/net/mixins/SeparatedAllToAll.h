#ifndef HEMELB_NET_MIXINS_SEPARATEDALLTOALL_H
#define HEMELB_NET_MIXINS_SEPARATEDALLTOALL_H

#include "net/mixins/StoringNet.h"

namespace hemelb{
  namespace net{

    class SeparatedAllToAll : public virtual StoringNet
    {
      void ReceiveAllToAll(){}
      void SendAllToAll(){}
      void WaitAllToAll();
    };
  }
}
#endif
