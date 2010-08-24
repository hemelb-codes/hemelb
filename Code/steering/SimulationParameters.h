#ifndef HEME_STEERING_SIMULATIONPARAMETERS_H
#define HEME_STEERING_SIMULATIONPARAMETERS_H
#ifndef NO_STEER

#include "io/XdrMemWriter.h"
#include "lb.h"

namespace heme
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
      heme::io::XdrMemWriter *paramWriter;

    };

  }
}

#endif // NO_STEER

#endif//__steering_simulationParameters_h_
