#ifndef HEMELB_UNITTESTS_NET_PHASED_MOCKITERATEDACTION_H
#define HEMELB_UNITTESTS_NET_PHASED_MOCKITERATEDACTION_H
#include <string>

#include "net/IteratedAction.h"
namespace hemelb
{
  namespace unittests
  {

    namespace net
    {
      using namespace hemelb::net;
      namespace phased
      {
        class MockIteratedAction : public IteratedAction
        {
          public:
            MockIteratedAction(const std::string & name):name(name){}
          private:
            std::string name;
        };
      }
    }
  }
}
#endif
