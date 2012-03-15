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

        LatticeTime GetTimeStep() const;
        LatticeTime Get0IndexedTimeStep() const;

        LatticeTime GetTimeStepsPassed() const;
        LatticeTime GetTotalTimeSteps() const;
        bool GetIsTerminating() const;
        bool GetDoRendering() const;
        Stability GetStability() const;

        PhysicalTime GetTime() const {return GetTimeStepLength()*Get0IndexedTimeStep();}
        PhysicalTime GetTimeStepLength() const {return TimeStepLength;}
        void DoubleTimeResolution();

        void Report(ctemplate::TemplateDictionary& dictionary);

      private:
        PhysicalTime TimeStepLength;
        LatticeTime TimeStep;
        LatticeTime TimeStepsGone;
        LatticeTime TotalTimeSteps;
        bool IsTerminating;
        bool DoRendering;
        Stability mStability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
