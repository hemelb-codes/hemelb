#include "constants.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    SimulationState::SimulationState(double timeStepLength, unsigned long totalTimeSteps) :
         TimeStepLength(timeStepLength), TimeStep(1), TimeStepsGone(1),
        TotalTimeSteps(totalTimeSteps), IsTerminating(false), DoRendering(false), mStability(Stable)
    {
    }

    void SimulationState::Increment()
    {
      ++TimeStepsGone;
      ++TimeStep;
    }

    void SimulationState::Reset()
    {
      TimeStep = 1;
      TimeStepsGone = 1;
    }

    void SimulationState::SetIsTerminating(bool value)
    {
      IsTerminating = value;
    }
    void SimulationState::SetDoRendering(bool value)
    {
      DoRendering = value;
    }
    void SimulationState::SetStability(Stability value)
    {
      mStability = value;
    }

    unsigned long SimulationState::GetTimeStep() const
    {
      return TimeStep;
    }

    unsigned long SimulationState::Get0IndexedTimeStep() const
    {
      return TimeStep - 1;
    }

    unsigned long SimulationState::GetTotalTimeSteps() const
    {
      return TotalTimeSteps;
    }

    unsigned long SimulationState::GetTimeStepsPassed() const
    {
      return TimeStepsGone;
    }

    bool SimulationState::GetIsTerminating() const
    {
      return IsTerminating;
    }
    bool SimulationState::GetDoRendering() const
    {
      return DoRendering;
    }
    Stability SimulationState::GetStability() const
    {
      return mStability;
    }

    void SimulationState::DoubleTimeResolution()
    {
      TotalTimeSteps *= 2;
      TimeStepLength /= 2.0;
    }

    void SimulationState::Report(ctemplate::TemplateDictionary& dictionary)
    {
      dictionary.SetFormattedValue("TIME_STEP_LENGTH","%lf",GetTimeStepLength());
      dictionary.SetIntValue("STEPS", GetTimeStepsPassed() - 1);
      dictionary.SetIntValue("TOTAL_TIME_STEPS", GetTotalTimeSteps());
    }
  }
}
