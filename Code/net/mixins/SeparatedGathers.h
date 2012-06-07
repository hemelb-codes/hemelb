#ifndef HEMELB_NET_MIXINS_SEPARATEDGATHERS_H
#define HEMELB_NET_MIXINS_SEPARATEDGATHERS_H
#include "net/BaseNet.h"
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
