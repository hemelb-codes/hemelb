#include "net/phased/StepManager.h"
namespace hemelb
{
  namespace net
  {
    namespace phased
    {

      StepManager::StepManager(Phase phases) :
          registry(phases)
      {
      }

      void StepManager::Register(Phase phase, steps::Step step, Concern & concern, MethodLabel method)
      {
        if (step==steps::BeginAll || step==steps::EndAll || step==steps::Reset) {
          phase=0; // special actions are always recorded in the phase zero registry
        }
        registry[phase].insert(BoundAction(step, Action(concern, method)));
      }

      void StepManager::RegisterIteratedActorSteps(Concern &concern, Phase phase)
      {
        for (int step = steps::BeginPhase; step <= steps::Reset; step++)
        {
          if (step == steps::Receive || step == steps::Send || step == steps::Wait)
          {
            continue;
          }
          // C++ makes it completely annoying to iterate over values in an enum
          Register(phase, static_cast<steps::Step>(step), concern, step);
        }
      }

      void StepManager::RegisterCommsSteps(Concern &concern, Phase phase)
      {
        Register(phase, steps::Send, concern, steps::Send);
        Register(phase, steps::Receive, concern, steps::Receive);
        Register(phase, steps::Wait, concern, steps::Wait);
      }

      unsigned int StepManager::ConcernCount() const
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

      unsigned int StepManager::ActionCount() const
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

      void StepManager::CallActionsForPhase(Phase phase)
      {
        for (int step = steps::BeginPhase; step <= steps::EndPhase; step++)
        {
          CallActionsForStep(static_cast<steps::Step>(step),phase);
        }
      }

      void StepManager::CallSpecialAction(steps::Step step)
      {
        // special actions are always recorded in the phase zero registry
        CallActionsForStep(static_cast<steps::Step>(step),0);
      }

      void StepManager::CallActions()
      {
        CallSpecialAction(steps::BeginAll);
        for (Phase phase=0;phase<registry.size();phase++){
          CallActionsForPhase(phase);
        }
        CallSpecialAction(steps::EndAll);
      }

      void StepManager::CallActionsForStep(steps::Step step,Phase phase)
      {

        std::pair<Registry::iterator, Registry::iterator> actionsForStep =
            registry[phase].equal_range(step);
        for (Registry::iterator action = actionsForStep.first; action != actionsForStep.second; ++action)
        {
          action->second.Call();
        }
      }
    }
  }
}
