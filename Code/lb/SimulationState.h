#ifndef SIMULATIONSTATE_H_
#define SIMULATIONSTATE_H_

#include <semaphore.h>

namespace hemelb
{
  namespace lb
  {
    struct SimulationState
    {
        int CycleId;
        int TimeStep;
        double IntraCycleTime;
        int IsTerminating;
        int DoRendering;
        int Stability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
