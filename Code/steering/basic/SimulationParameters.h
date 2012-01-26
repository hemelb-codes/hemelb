#ifndef HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H
#define HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H

#include "io/writers/xdr/XdrMemWriter.h"
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
        void pack(io::writers::xdr::XdrWriter* writer);

      private:

    };

  }
}

#endif // HEMELB_STEERING_BASIC_SIMULATIONPARAMETERS_H
