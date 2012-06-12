#ifndef HEMELB_UNITTESTS_NET_PHASED_MOCKITERATEDACTION_H
#define HEMELB_UNITTESTS_NET_PHASED_MOCKITERATEDACTION_H
#include <string>
#include <sstream>

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
            MockIteratedAction(const std::string & name) :
                name(name), calls()
            {
            }
            void PreSend()
            {
              calls << "PreSend, " << std::flush;
            }
            void PreReceive()
            {
              calls << "PreReceive, " << std::flush;
            }
            void PostReceive()
            {
              calls << "PostReceive, " << std::flush;
            }
            void EndIteration()
            {
              calls << "EndIteration, " << std::flush;
            }
            void RequestComms()
            {
              calls << "RequestComms, " << std::flush;
            }
            std::string CallsSoFar(){
              return calls.str();
            }

          private:
            std::string name;
            std::stringstream calls;
        };
      }
    }
  }
}
#endif
