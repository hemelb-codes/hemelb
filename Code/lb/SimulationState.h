#ifndef SIMULATIONSTATE_H_
#define SIMULATIONSTATE_H_

#include <semaphore.h>

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

    struct SimulationState
    {
      public:
        SimulationState(unsigned long StepsPerCycle, unsigned long numCycles);

        void Increment();
        void Reset();
        void SetIsTerminating(bool value);
        void SetDoRendering(bool value);
        void SetStability(Stability value);

        unsigned long GetCycleId() const;
        unsigned long GetTimeStep() const;
        unsigned long Get0IndexedTimeStep() const;
        unsigned long GetTimeStepsPerCycle() const;
        unsigned long GetNumberOfCycles() const;
        unsigned long GetTimeStepsPassed() const;
        unsigned long GetTotalTimeSteps() const;
        double GetIntraCycleTime() const;
        bool GetIsTerminating() const;
        bool GetDoRendering() const;
        Stability GetStability() const;

        void DoubleTimeResolution();

      private:
        unsigned long CycleId;
        unsigned long TimeStep;
        unsigned long TimeStepsGone;
        unsigned long TotalTimeSteps;
        unsigned long TimeStepsPerCycle;
        unsigned long NumberOfCycles;
        bool IsTerminating;
        bool DoRendering;
        Stability mStability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
