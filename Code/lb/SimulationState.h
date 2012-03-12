#ifndef HEMELB_LB_SIMULATIONSTATE_H
#define HEMELB_LB_SIMULATIONSTATE_H

#include "reporting/Reportable.h"

namespace hemelb
{
  namespace lb
  {
    enum Stability
    {
      Unstable = 0,
      Stable = 1,
      StableAndConverged = 2
    };

    class SimulationState : public reporting::Reportable
    {
      public:
        SimulationState(double timeStepLength, unsigned long totalTimeSteps);

        void Increment();
        void Reset();
        void SetIsTerminating(bool value);
        void SetDoRendering(bool value);
        void SetStability(Stability value);

        unsigned long GetTimeStep() const;
        unsigned long Get0IndexedTimeStep() const;

        unsigned long GetTimeStepsPassed() const;
        unsigned long GetTotalTimeSteps() const;
        bool GetIsTerminating() const;
        bool GetDoRendering() const;
        Stability GetStability() const;

        double GetTimeStepLength() const {return TimeStepLength;}
        void DoubleTimeResolution();

        void Report(ctemplate::TemplateDictionary& dictionary);

      private:
        double TimeStepLength;
        unsigned long TimeStep;
        unsigned long TimeStepsGone;
        unsigned long TotalTimeSteps;
        bool IsTerminating;
        bool DoRendering;
        Stability mStability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
