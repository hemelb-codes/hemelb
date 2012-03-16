#include "constants.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace lb
  {

    SimulationState::SimulationState(double atimeStepLength, unsigned long atotalTimeSteps) :
        timeStepLength(atimeStepLength), timeStep(1), totalTimeSteps(atotalTimeSteps), isTerminating(false), isRendering(false), stability(Stable)
    {
    }

    void SimulationState::Increment()
    {
      ++timeStep;
    }

    void SimulationState::Reset()
    {
      timeStep = 1;
    }

    void SimulationState::SetIsTerminating(bool value)
    {
      isTerminating = value;
    }
    void SimulationState::SetIsRendering(bool value)
    {
      isRendering = value;
    }
    void SimulationState::SetStability(Stability value)
    {
      stability = value;
    }

    unsigned long SimulationState::GetTimeStep() const
    {
      return timeStep;
    }

    unsigned long SimulationState::Get0IndexedTimeStep() const
    {
      return timeStep - 1;
    }

    unsigned long SimulationState::GetTotalTimeSteps() const
    {
      return totalTimeSteps;
    }

    bool SimulationState::IsTerminating() const
    {
      return isTerminating;
    }
    bool SimulationState::IsRendering() const
    {
      return isRendering;
    }
    Stability SimulationState::GetStability() const
    {
      return stability;
    }

    void SimulationState::DoubleTimeResolution()
    {
      totalTimeSteps *= 2;
      timeStepLength /= 2.0;
    }

    void SimulationState::Report(ctemplate::TemplateDictionary& dictionary)
    {
      dictionary.SetFormattedValue("TIME_STEP_LENGTH", "%lf", GetTimeStepLength());
      dictionary.SetIntValue("STEPS", GetTimeStep() - 1);
      dictionary.SetIntValue("TOTAL_TIME_STEPS", GetTotalTimeSteps());
      if (stability == lb::Unstable)
      {
        dictionary.AddSectionDictionary("UNSTABLE");
      }
    }
  }
}
