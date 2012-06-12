#include "net/IteratedAction.h"
#include "net/phased/StepManager.h"

namespace hemelb
{
  namespace net
  {
    IteratedAction::~IteratedAction()
    {

    }

    bool IteratedAction::CallAction(int action)
    {
      switch (static_cast<phased::steps::Step>(action))
      {
        case phased::steps::BeginPhase:
          RequestComms();
          return true;
        case phased::steps::PreSend:
          PreSend();
          return true;
        case phased::steps::PreWait:
          PreReceive();
          return true;
        case phased::steps::EndPhase:
          PostReceive();
          return true;
        case phased::steps::EndAll:
          EndIteration();
          return true;
        default:
          return false;
      }
    }

    void IteratedAction::RequestComms()
    {

    }

    void IteratedAction::PreSend()
    {

    }

    void IteratedAction::PreReceive()
    {

    }

    void IteratedAction::PostReceive()
    {

    }

    void IteratedAction::EndIteration()
    {

    }

    void IteratedAction::Reset()
    {

    }
  }
}
