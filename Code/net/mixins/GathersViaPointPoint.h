#ifndef HEMELB_NET_MIXINS_GATHERSVIAPOINTPOINT_H
#define HEMELB_NET_MIXINS_GATHERSVIAPOINTPOINT_H
#include "net/BaseNet.h"
#include "net/mixins/StoringNet.h"
namespace hemelb{
  namespace net{
    class GathersViaPointPoint : public virtual StoringNet
    {
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
