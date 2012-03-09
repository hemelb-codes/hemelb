#include "constants.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    SimulationState::SimulationState(unsigned long StepsPerCycle, unsigned long numCycles)
    {
      CycleId = 1;
      TimeStep = 1;
      TimeStepsGone = 1;
      TimeStepsPerCycle = StepsPerCycle;
      NumberOfCycles = numCycles;
      TotalTimeSteps = numCycles * StepsPerCycle;
      IsTerminating = false;
      DoRendering = false;
      mStability = Stable;
    }

    void SimulationState::Increment()
    {
      ++TimeStepsGone;

      if (TimeStep >= TimeStepsPerCycle)
      {
        ++CycleId;
        TimeStep = 1;
      }
      else
      {
        ++TimeStep;
      }
    }

    void SimulationState::Reset()
    {
      CycleId = 1;
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

    unsigned long SimulationState::GetCycleId() const
    {
      return CycleId;
    }

    unsigned long SimulationState::GetTimeStep() const
    {
      return TimeStep;
    }

    unsigned long SimulationState::Get0IndexedTimeStep() const
    {
      return TimeStep - 1;
    }

    unsigned long SimulationState::GetTimeStepsPerCycle() const
    {
      return TimeStepsPerCycle;
    }

    unsigned long SimulationState::GetNumberOfCycles() const
    {
      return NumberOfCycles;
    }

    unsigned long SimulationState::GetTotalTimeSteps() const
    {
      return TotalTimeSteps;
    }

    unsigned long SimulationState::GetTimeStepsPassed() const
    {
      return TimeStepsGone;
    }

    double SimulationState::GetIntraCycleTime() const
    {
      return hemelb::PULSATILE_PERIOD_s * (double) TimeStep / (double) TimeStepsPerCycle;
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
      TimeStepsPerCycle *= 2;
    }

    void SimulationState::Report(ctemplate::TemplateDictionary& dictionary)
    {
      // Note that CycleId is 1-indexed and will have just been incremented when we finish.
      dictionary.SetIntValue("CYCLES", GetCycleId() - 1);
      dictionary.SetIntValue("STEPS", GetTimeStepsPassed() - 1);
      dictionary.SetIntValue("STEPS_PER_CYCLE", GetTimeStepsPerCycle());
    }
  }
}
