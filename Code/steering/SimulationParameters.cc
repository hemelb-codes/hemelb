#ifndef NO_STEER

#include "steering/SimulationParameters.h"
#include "vis/rt.h"
#include "io/XdrMemWriter.h"

using namespace std;

extern int cycle_id;
extern int time_step;
extern double intra_cycle_time;

heme::steering::SimulationParameters::SimulationParameters() {
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
  paramWriter = new heme::io::XdrMemWriter(params, paramsSizeB);

}

void heme::steering::SimulationParameters :: collectGlobalVals() {
  
  this->pressureMin = heme::vis::pressure_min;
  this->pressureMax = heme::vis::pressure_max;
  this->velocityMin = heme::vis::velocity_min;
  this->velocityMax = heme::vis::velocity_max;
  this->stressMax = heme::vis::stress_max;
  this->timeStep = time_step;
  this->time = intra_cycle_time;
  this->cycle = cycle_id;
  this->nInlets = heme::vis::inlets;
  
  this->mousePressure = heme::vis::mouse_pressure;
  this->mouseStress = heme::vis::mouse_stress;
  
}

heme::steering::SimulationParameters::~SimulationParameters() {
  delete paramWriter;
  delete[] inletAvgVel;
  // TODO: find out if there's a good reason this isn't deleted
  // delete[] params;
}

char* heme::steering::SimulationParameters::pack() {
  heme::io::XdrMemWriter& paramWriter = *(this->paramWriter);
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

  heme::vis::mouse_pressure = -1.0;
  heme::vis::mouse_stress = -1.0;

  return params;
}

#endif // NO_STEER
