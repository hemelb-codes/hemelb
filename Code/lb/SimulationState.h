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
        sem_t Rendering;

        SimulationState()
        {
          sem_init(&Rendering, 0, 1);
        }

        ~SimulationState()
        {
          sem_destroy(&Rendering);
        }
    };
  }
}

#endif /* SIMULATIONSTATE_H_ */
