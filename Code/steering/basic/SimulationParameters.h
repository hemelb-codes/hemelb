#ifndef HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H
#define HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H

#include "io/XdrMemWriter.h"
#include "lb.h"

namespace hemelb
{
  namespace steering
  {

    class SimulationParameters {
  
    public:
  
      double pressureMin;
      double pressureMax;
      double velocityMin;
      double velocityMax;
      double stressMin;
      double stressMax;
      int timeStep;
      double time;
      int cycle;
      int nInlets;
      double* inletAvgVel;
      double mousePressure;
      double mouseStress;

      char* params;
      u_int paramsSizeB;

      SimulationParameters();
      ~SimulationParameters();
      char* pack();
      void collectGlobalVals(LBM* lbm);

    private:
      io::XdrMemWriter *paramWriter;

    };

  }
}


#endif // HEMELB_STEERING_ON_SIMULATIONPARAMETERS_H
