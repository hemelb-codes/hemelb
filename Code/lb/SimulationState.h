#ifndef SIMULATIONSTATE_H_
#define SIMULATIONSTATE_H_

#include <semaphore.h>

namespace hemelb
{
  namespace lb
  {
    struct SimulationState
    {
        unsigned long CycleId;
        unsigned long TimeStep;
        unsigned long TimeStepsPerCycle;
        double IntraCycleTime;
        int IsTerminating;
        int DoRendering;
        int Stability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
