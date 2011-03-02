#ifndef SIMULATIONSTATE_H_
#define SIMULATIONSTATE_H_

#include <semaphore.h>

namespace hemelb
{
  namespace lb
  {
    struct SimulationState
    {
        long CycleId;
        unsigned int TimeStep;
        double IntraCycleTime;
        int IsTerminating;
        int DoRendering;
        int Stability;
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
