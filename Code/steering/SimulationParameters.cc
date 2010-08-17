#ifndef NO_STEER

#include "steering/SimulationParameters.h"
#include "vis/rt.h"
#include "io/xdrMemWriter.h"

using namespace std;

extern int cycle_id;
extern int time_step;
extern double intra_cycle_time;

SimulationParameters::SimulationParameters() {
  // C'tor initialises to the following defaults.
  
  pressureMin = 0.001;
  pressureMax = 1.0;
  velocityMin = 0.0;
  velocityMax = 1.0;
  stressMin = 0.0;
  stressMax = 1.0;
  timeStep = 0;
  time = 0.0;
  cycle = 0;
  nInlets = 3;
  mousePressure = -1.0;
  mouseStress = -1.0;

  inletAvgVel = new double[nInlets];

  for(int i=0; i<nInlets; i++)
    inletAvgVel[i] = 1.0;

  // Assumption here is that sizeof(char) is 1B;
  paramsSizeB = 3 * sizeof(int);
  paramsSizeB += 9 * sizeof(double);
  // params_bytes += n_inlets * sizeof(double);

  params = new char[paramsSizeB];
  paramWriter = new io::XdrMemWriter(params, paramsSizeB);

}

void SimulationParameters :: collectGlobalVals() {
  
  this->pressureMin = vis::pressure_min;
  this->pressureMax = vis::pressure_max;
  this->velocityMin = vis::velocity_min;
  this->velocityMax = vis::velocity_max;
  this->stressMax = vis::stress_max;
  this->timeStep = time_step;
  this->time = intra_cycle_time;
  this->cycle = cycle_id;
  this->nInlets = vis::inlets;
  
  this->mousePressure = vis::mouse_pressure;
  this->mouseStress = vis::mouse_stress;
  
}

SimulationParameters::~SimulationParameters() {
  delete paramWriter;
  delete[] inletAvgVel;
  // TODO: find out if there's a good reason this isn't deleted
  // delete[] params;
}

char* SimulationParameters::pack() {
  io::XdrMemWriter& paramWriter = *(this->paramWriter);
  paramWriter << pressureMin;
  paramWriter << pressureMax;

  paramWriter << velocityMin;
  paramWriter << velocityMax;

  paramWriter << stressMin;
  paramWriter << stressMax;

  paramWriter << timeStep;

  paramWriter << time;

  paramWriter << cycle;
  paramWriter << nInlets;
  
  //	for(int i=0; i<n_inlets; i++) xdr_double(&xdr_params, &inlet_avg_vel[i]);

  paramWriter << mousePressure;
  paramWriter << mouseStress;

  vis::mouse_pressure = -1.0;
  vis::mouse_stress = -1.0;

  return params;
}

#endif // NO_STEER
