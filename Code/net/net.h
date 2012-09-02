#ifndef HEMELB_NET_NET_H
#define HEMELB_NET_NET_H

#include "net/BaseNet.h"
#include "net/mixins/mixins.h"
namespace hemelb
{
  namespace net
  {
    class Net : public SeparatedPointPoint,
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

