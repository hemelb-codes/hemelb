#ifndef __steering_simulationParameters_h_
#define __steering_simulationParameters_h_
#ifndef NO_STEER

#include "io/xdrMemWriter.h"

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

  XDR xdr_sim_params;
  char* params;
  u_int paramsSizeB;

  SimulationParameters();
  ~SimulationParameters();
  char* pack();
  void collectGlobalVals();

 private:
  io::XdrMemWriter *paramWriter;

};

#endif // NO_STEER

#endif//__steering_simulationParameters_h_
