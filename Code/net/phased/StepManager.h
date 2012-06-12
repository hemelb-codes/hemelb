#ifndef HEMELB_NET_PHASED_STEPMANAGER_H
#define HEMELB_NET_PHASED_STEPMANAGER_H
#include <vector>
#include <map>
#include <set>
#include "net/IteratedAction.h"
namespace hemelb
{
  namespace net
  {
    typedef IteratedAction Concern; // temporary hack until we implement concern.
    namespace phased
    {
      enum Step
       {
         RequestComms,
         Receive,
         PreSend,
         Send,
         PreReceive,
         Wait,
         PostReceive
       } ;

      class StepManager
      {
        public:


          typedef unsigned int Phase;
          typedef unsigned int MethodLabel;

          class Action
          {
            public:
              Concern & concern;
              MethodLabel method;
              Action(Concern &concern,MethodLabel method):concern(concern),method(method){}
          };

          typedef std::pair<Step, Action> BoundAction;
          typedef std::multimap<Step, Action> Registry;

          StepManager(Phase phases=1) :
              registry(phases)
          {
          }
          void Register(Phase phase, Step step, Concern & concern, MethodLabel method)
          {
            registry[phase].insert(BoundAction(step, Action(concern, method)));
          }
          void RegisterAllSteps(Concern &concern, Phase phase = 0)
          {
            for (int step = RequestComms; step <= PostReceive; step++)
            {
              // C++ makes it completely annoying to iterate over values in an enum
              Register(phase, static_cast<Step>(step), concern, step);
            }
          }

          unsigned int ConcernCount() const
          {
            std::set<Concern *> concerns;
            for (std::vector<Registry>::const_iterator phase = registry.begin(); phase < registry.end(); phase++)
            {
              for (Registry::const_iterator boundAction = phase->begin(); boundAction != phase->end(); boundAction++)
              {
                concerns.insert(&boundAction->second.concern);
              }
            }
            return concerns.size();
          }

          unsigned int ActionCount() const
          {
            unsigned int total = 0;
            for (std::vector<Registry>::const_iterator phase = registry.begin(); phase < registry.end(); phase++)
            {
              for (Registry::const_iterator boundAction = phase->begin(); boundAction != phase->end(); boundAction++)
              {
                total++;
              }
            }
            return total;
          }
        private:
          std::vector<Registry> registry;

      };
    }
  }
}
#endif //ONCE
