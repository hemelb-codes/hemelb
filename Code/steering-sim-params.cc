#ifndef NO_STEER

#include <stdio.h>
#include <stdlib.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "steering-sim-params.h"
#include "vis/rt.h"

using namespace std;

extern int cycle_id;
extern int time_step;
extern double intra_cycle_time;

simulationParameters :: simulationParameters() {

	sim_pressure_min = 0.001;
	sim_pressure_max = 1.0;
	sim_velocity_min = 0.0;
	sim_velocity_max = 1.0;
	sim_stress_min = 0.0;
	sim_stress_max = 1.0;
	sim_time_step = 0;
	sim_time = 0.0;
	sim_cycle = 0;
	sim_n_inlets = 3;
	sim_mouse_pressure = -1.0;
	sim_mouse_stress = -1.0;

	sim_inlet_avg_vel = new double[sim_n_inlets];

	for(int i=0; i<sim_n_inlets; i++) sim_inlet_avg_vel[i] = 1.0;

	// Assumption here is that sizeof(char) is 1B;
	sim_params_bytes = 3 * sizeof(int);
	sim_params_bytes += 9 * sizeof(double);
//	sim_params_bytes += sim_n_inlets * sizeof(double);

	sim_params = new char[sim_params_bytes];

	xdrmem_create(&xdr_sim_params, sim_params, sim_params_bytes, XDR_ENCODE);

}

void simulationParameters :: collectGlobalVals() {
  
  this->sim_pressure_min = vis::pressure_min;
  this->sim_pressure_max = vis::pressure_max;
  this->sim_velocity_min = vis::velocity_min;
  this->sim_velocity_max = vis::velocity_max;
  this->sim_stress_max = vis::stress_max;
  this->sim_time_step = time_step;
  this->sim_time = intra_cycle_time;
  this->sim_cycle = cycle_id;
  this->sim_n_inlets = vis::inlets;
  
  //for(int i=0; i<sim_n_inlets; i++) this->sim_inlet_avg_vel[i] = lbm_inlet_flux[i];
  // for(int i=0; i<sim_n_inlets; i++) printf("avg vel %0.3f\n", this->sim_inlet_avg_vel[i]);
  
  this->sim_mouse_pressure = vis::mouse_pressure;
  this->sim_mouse_stress = vis::mouse_stress;

}

simulationParameters :: ~simulationParameters() {
  xdr_destroy(&xdr_sim_params);
  delete[] sim_inlet_avg_vel;
  //	delete[] sim_params;
}

char* simulationParameters :: pack() {

  xdr_double(&xdr_sim_params, &sim_pressure_min);
  xdr_double(&xdr_sim_params, &sim_pressure_max);

  xdr_double(&xdr_sim_params, &sim_velocity_min);
  xdr_double(&xdr_sim_params, &sim_velocity_max);

  xdr_double(&xdr_sim_params, &sim_stress_min);
  xdr_double(&xdr_sim_params, &sim_stress_max);

  xdr_int(&xdr_sim_params, &sim_time_step);

  xdr_double(&xdr_sim_params, &sim_time);

  xdr_int(&xdr_sim_params, &sim_cycle);
  xdr_int(&xdr_sim_params, &sim_n_inlets);

  //	for(int i=0; i<sim_n_inlets; i++) xdr_double(&xdr_sim_params, &sim_inlet_avg_vel[i]);

  xdr_double(&xdr_sim_params, &sim_mouse_pressure);
  xdr_double(&xdr_sim_params, &sim_mouse_stress);

  vis::mouse_pressure = -1.0;
  vis::mouse_stress = -1.0;

  return sim_params;
}

#endif // NO_STEER
