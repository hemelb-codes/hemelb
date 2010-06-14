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

  // TODO Ideally this block wouldn't exist.
  sim_pressure_min = 0.0;
  sim_pressure_max = 0.0;
  sim_velocity_min = 0.0;
  sim_velocity_max = 0.0;
  sim_stress_min = 0.0;
  sim_stress_max = 0.0;
  sim_time_step = 0;
  sim_time = 0.0;
  sim_cycle = 0;
  sim_n_inlets = 0;
  sim_mouse_pressure = -1.0;
  sim_mouse_stress = -1.0;


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
  sim_time_step = time_step;
  sim_time = intra_cycle_time;
  sim_cycle = cycle_id;
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

  sim_mouse_pressure = -1.0;
  sim_mouse_stress = -1.0;

  return sim_params;
}

// Number of bytes required to pack the simulation parameters
u_int simulationParameters :: getPackedSizeInBytes()
{
  // Assuming that sizeof(char) is 1B;
  return 3 * sizeof(int) + 9 * sizeof(double);
}

// Setter methods
void simulationParameters :: set_Min_Sim_Pressure(double new_min_pressure) {sim_pressure_min = new_min_pressure;}
void simulationParameters :: set_Max_Sim_Pressure(double new_max_pressure) {sim_pressure_max = new_max_pressure;}
void simulationParameters :: set_Min_Sim_Velocity(double new_min_velocity) {sim_velocity_min = new_min_velocity;}
void simulationParameters :: set_Max_Sim_Velocity(double new_max_velocity) {sim_velocity_max = new_max_velocity;}
void simulationParameters :: set_Min_Sim_Stress(double new_min_stress) {sim_stress_min = new_min_stress;}
void simulationParameters :: set_Max_Sim_Stress(double new_max_stress) {sim_stress_max = new_max_stress;}

void simulationParameters :: set_Sim_Inlets(int n_inlets) {sim_n_inlets = n_inlets;}

void simulationParameters :: set_Sim_Mouse_Pressure(double new_mouse_pressure) {sim_mouse_pressure = new_mouse_pressure;}
void simulationParameters :: set_Sim_Mouse_Stress(double new_mouse_stress) {sim_mouse_stress = new_mouse_stress;}

// Accessor methods
double simulationParameters :: get_Min_Sim_Pressure() {return sim_pressure_min;}
double simulationParameters :: get_Max_Sim_Pressure() {return sim_pressure_max;}
double simulationParameters :: get_Min_Sim_Velocity() {return sim_velocity_min;}
double simulationParameters :: get_Max_Sim_Velocity() {return sim_velocity_max;}
double simulationParameters :: get_Min_Sim_Stress() {return sim_stress_min;}
double simulationParameters :: get_Max_Sim_Stress() {return sim_stress_max;}

#endif // NO_STEER
