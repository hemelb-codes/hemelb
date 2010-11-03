#include "io/XdrMemWriter.h"

#include "vis/Control.h"
#include "steering/basic/SimulationParameters.h"

using namespace std;

hemelb::steering::SimulationParameters::SimulationParameters()
{
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

  for (int i = 0; i < nInlets; i++)
    inletAvgVel[i] = 1.0;

  // Assumption here is that sizeof(char) is 1B;
  paramsSizeB = 3 * sizeof(int);
  paramsSizeB += 9 * sizeof(double);
  // params_bytes += n_inlets * sizeof(double);

  params = new char[paramsSizeB];
  paramWriter = new io::XdrMemWriter(params, paramsSizeB);

}

void hemelb::steering::SimulationParameters::collectGlobalVals(LBM* lbm,
                                                               lb::SimulationState *iSimState)
{
  this->pressureMin = lbm->GetMinPhysicalPressure();
  this->pressureMax = lbm->GetMaxPhysicalPressure();
  this->velocityMin = lbm->GetMinPhysicalVelocity();
  this->velocityMax = lbm->GetMaxPhysicalVelocity();
  this->stressMax = lbm->GetMaxPhysicalStress();
  this->timeStep = iSimState->TimeStep;
  this->time = iSimState->IntraCycleTime;
  this->cycle = iSimState->CycleId;
  this->nInlets = lbm->inlets;

  this->mousePressure = vis::controller->mouse_pressure;
  this->mouseStress = vis::controller->mouse_stress;

}

hemelb::steering::SimulationParameters::~SimulationParameters()
{
  delete paramWriter;
  delete[] inletAvgVel;
  // TODO: find out if there's a good reason this isn't deleted
  // delete[] params;
}

char* hemelb::steering::SimulationParameters::pack()
{
  io::XdrMemWriter& paramWriter = * (this->paramWriter);
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

  vis::controller->mouse_pressure = -1.0;
  vis::controller->mouse_stress = -1.0;

  return params;
}

