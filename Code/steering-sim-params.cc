#ifndef NO_STEER

#include <stdio.h>
#include <stdlib.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "config.h"
#include "steering-sim-params.h"

using namespace std;

extern int cycle_id;
extern int time_step;
extern double intra_cycle_time;

simulationParameters :: simulationParameters() 
{
  int sim_params_length = getPackedSizeInBytes();
  sim_params = new char[sim_params_length];
  xdrmem_create(&xdr_sim_params, sim_params, sim_params_length, XDR_ENCODE);
}

simulationParameters :: ~simulationParameters() 
{
  xdr_destroy(&xdr_sim_params);
  delete[] sim_params;
}

void simulationParameters :: collectGlobalVals() 
{
  sim_pressure_min = vis_pressure_min;
  sim_pressure_max = vis_pressure_max;
  sim_velocity_min = vis_velocity_min;
  sim_velocity_max = vis_velocity_max;
  sim_stress_max = vis_stress_max;
  sim_time_step = time_step;
  sim_time = intra_cycle_time;
  sim_cycle = cycle_id;
  sim_n_inlets = vis_inlets;
  sim_mouse_pressure = vis_mouse_pressure;
  sim_mouse_stress = vis_mouse_stress;
}

char* simulationParameters :: pack() 
{
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

  xdr_double(&xdr_sim_params, &sim_mouse_pressure);
  xdr_double(&xdr_sim_params, &sim_mouse_stress);

  vis_mouse_pressure = -1.0;
  vis_mouse_stress = -1.0;

  return sim_params;
}

// Number of bytes required to pack the simulation parameters
u_int simulationParameters :: getPackedSizeInBytes()
{
  // Assuming that sizeof(char) is 1B;
  return 3 * sizeof(int) + 9 * sizeof(double);
}

#endif // NO_STEER
