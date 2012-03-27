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
        void SetIsRendering(bool value);
        void SetStability(Stability value);

        LatticeTime GetTimeStep() const;
        LatticeTime Get0IndexedTimeStep() const;
        LatticeTime GetTotalTimeSteps() const;
        bool IsTerminating() const;
        bool IsRendering() const;
        Stability GetStability() const;

        PhysicalTime GetTime() const {return GetTimeStepLength()*Get0IndexedTimeStep();}
        PhysicalTime GetTimeStepLength() const {return timeStepLength;}
        void DoubleTimeResolution();

        void Report(ctemplate::TemplateDictionary& dictionary);

      private:
        PhysicalTime timeStepLength;
        LatticeTime timeStep;
        LatticeTime totalTimeSteps;
        bool isTerminating;
        bool isRendering;
        Stability stability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
