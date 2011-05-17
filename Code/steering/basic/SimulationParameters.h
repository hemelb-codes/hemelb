#ifndef HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H
#define HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H

#include "io/XdrMemWriter.h"
#include "lb/SimulationState.h"

namespace hemelb
{
  namespace steering
  {

    class SimulationParameters
    {

      public:

        int timeStep;
        double time;
        int cycle;
        int nInlets;
        double mousePressure;
        double mouseStress;

        static const u_int paramsSizeB = 3 * sizeof(int) + 3 * sizeof(double);

        SimulationParameters();
        ~SimulationParameters();
        void pack(io::XdrWriter* writer);

      private:

    };

  }
}

#endif // HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H
