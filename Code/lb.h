#ifndef __lb_h_
#define __lb_h_

#include "net.h"

struct LBM
{
  char *system_file_name;
  
  double tau, viscosity;
  double voxel_size;
  double omega;
  
  int total_fluid_sites;
  int site_min_x, site_min_y, site_min_z;
  int site_max_x, site_max_y, site_max_z;
  int inlets, outlets;
  int cycles_max;
  int period;
  int conv_freq;
  
  float *block_density;
  
  int *block_map;
};

double lbmConvertPressureToLatticeUnits (double pressure, LBM *lbm);
double lbmConvertPressureToPhysicalUnits (double pressure, LBM *lbm);
double lbmConvertPressureGradToLatticeUnits (double pressure_grad, LBM *lbm);
double lbmConvertPressureGradToPhysicalUnits (double pressure_grad, LBM *lbm);
double lbmConvertVelocityToLatticeUnits (double velocity, LBM *lbm);
double lbmConvertVelocityToPhysicalUnits (double velocity, LBM *lbm);
double lbmConvertStressToLatticeUnits (double stress, LBM *lbm);
double lbmConvertStressToPhysicalUnits (double stress, LBM *lbm);
void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[]);
void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[]);
void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z);
void lbmStress (double f[], double *stress);
void lbmStress (double density, double f[], double nor[], double *stress);
void lbmInitMinMaxValues (void);
void lbmUpdateMinMaxValues (double density, double velocity, double stress);
void lbmCalculateBC (double f[], unsigned int site_data, double *density, double *vx, double *vy, double *vz, double f_neq[]);
int lbmCollisionType (unsigned int site_data);
void lbmInit (char *system_file_name, LBM *lbm, Net *net);
void lbmSetInitialConditions (LBM *lbm, Net *net);
void lbmUpdateFlowField (int perform_rt, int i, double density, double vx, double vy, double vz, double f_neq[]);
void lbmUpdateFlowFieldConv (int perform_rt, int i, double density, double vx, double vy, double vz, double f_neq[]);
int lbmCycle (int perform_rt, LBM *lbm, Net *net);
int lbmCycle (int cycle_id, int time_step, int perform_rt, LBM *lbm, Net *net);
void lbmCalculateFlowFieldValues (LBM *lbm);
int lbmIsUnstable (Net *net);
void lbmRestart (LBM *lbm, Net *net);
void lbmEnd (void);
void lbmReadConfig (LBM *lbm, Net *net);
double lbmCalculateTau (LBM *lbm);
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net);

void lbmWriteConfig (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmWriteConfigASCII (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmUpdateBoundaryDensities (int cycle_id, int time_step, LBM *lbm);
void lbmUpdateInletVelocities (int time_step, LBM *lbm, Net *net);



extern double lbm_stress_type;
extern double lbm_stress_par;
extern double lbm_density_min, lbm_density_max;
extern double lbm_velocity_min, lbm_velocity_max;
extern double lbm_stress_min, lbm_stress_max;
extern double *lbm_average_inlet_velocity;
extern double *lbm_peak_inlet_velocity;
extern double *lbm_inlet_normal;
extern long int *lbm_inlet_count;

extern int lbm_terminate_simulation;



extern double *inlet_density;
extern double *inlet_density_avg, *inlet_density_amp, *inlet_density_phs;
extern double *outlet_density;
extern double *outlet_density_avg, *outlet_density_amp, *outlet_density_phs;

extern int is_inlet_normal_available;


// TODO Judging by the name, these shouldn't be in here.

int netFindTopology (Net *net, int *depths);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

// TODO: prob don't belong here... 3 variables needed for convergence-enabled simulations
extern double conv_error;
extern int cycle_tag, check_conv;
#endif //__lb_h_
