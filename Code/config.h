// In this file all the declarations and some simple functions are reported.

#ifndef __config_h__
#define __config_h__

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

extern int STABLE;
extern int UNSTABLE;

// the constants needed to define the configuration of the lattice
// sites follow

extern unsigned int SOLID_TYPE;
extern unsigned int FLUID_TYPE;
extern unsigned int INLET_TYPE;
extern unsigned int OUTLET_TYPE;
extern unsigned int NULL_TYPE;

extern unsigned int BOUNDARIES;
extern unsigned int INLET_BOUNDARY;
extern unsigned int OUTLET_BOUNDARY;
extern unsigned int WALL_BOUNDARY;
extern unsigned int CHARACTERISTIC_BOUNDARY;

extern unsigned int SITE_TYPE_BITS;
extern unsigned int BOUNDARY_CONFIG_BITS;
extern unsigned int BOUNDARY_DIR_BITS;
extern unsigned int BOUNDARY_ID_BITS;

extern unsigned int BOUNDARY_CONFIG_SHIFT ;   // SITE_TYPE_BITS;
extern unsigned int BOUNDARY_DIR_SHIFT;  // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
extern unsigned int BOUNDARY_ID_SHIFT;  // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

extern unsigned int SITE_TYPE_MASK;         // ((1U << SITE_TYPE_BITS) - 1U);
extern unsigned int BOUNDARY_CONFIG_MASK;   // ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
extern unsigned int BOUNDARY_DIR_MASK;  //((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
extern unsigned int BOUNDARY_ID_MASK;  // ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT
extern unsigned int CHARACTERISTIC_MASK;

// square of the speed of sound

extern double Cs2;


// parameters related to the lattice directions

extern int e_x[];
extern int e_y[];
extern int e_z[];
extern int inv_dir[];


#define MACROSCOPIC_PARS   5

#define NEIGHBOUR_PROCS_MAX        52
#define SHARED_DISTRIBUTIONS_MAX   4000000


// data structures useful to define the simulation set-up, construct
// the system and partition it

struct BlockLocation
{
  short int i, j, k;
};


struct LBM
{
  char *system_file_name;
  char *parameters_file_name;
  char *checkpoint_file_name;
  
  double tau;
  double viscosity;
  double omega, stress_par;
  double lattice_to_system;
  double *inlet_density;
  double *outlet_density;
  double tolerance;
  
  int sites_x, sites_y, sites_z;
  int blocks_x, blocks_y, blocks_z;
  int blocks;
  int block_size;
  int shift;
  int sites_in_a_block;
  int total_fluid_sites;
  int site_min_x, site_min_y, site_min_z;
  int site_max_x, site_max_y, site_max_z;
  int inlets, outlets;
  int time_steps_max;
  int checkpoint_frequency, convergence_frequency;
  int is_checkpoint;
  
  short int *fluid_sites_per_block;
};


struct NeighProc
{
  int id;
  int fs;
  
  short int *f_data;
  
  double *f_to_send;
  double *f_to_recv;
  
  int *f_send_id;
  int *f_recv_iv;
};


struct DataBlock
{
  unsigned int *site_data;
};


struct Velocity
{
  double x, y, z;
};


struct Net
{
  int sites_x, sites_y, sites_z;
  int blocks_x, blocks_y, blocks_z;
  int blocks;
  int block_size;
  int shift;
  int sites_in_a_block;
  int id;
  int procs, machines;
  int inter_m_neigh_procs, neigh_procs;
  int err;
  int my_inter_sites, my_inner_sites;
  
  int *machine_id;
  int *procs_per_machine;
  
  short int *proc_id;
  
  unsigned int *site_data;
  
  DataBlock *data_block;
  DataBlock *map_block;
  
  NeighProc inter_m_neigh_proc[NEIGHBOUR_PROCS_MAX];
  NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
  MPI_Status status[4];
  
  MPI_Request **req;
  
  double timing[7];
};


extern double *f_old, *f_new;

extern int *f_id;

extern Velocity *vel;


extern short int f_data[4*SHARED_DISTRIBUTIONS_MAX];

extern double f_to_send[SHARED_DISTRIBUTIONS_MAX];
extern double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

extern int f_send_id[SHARED_DISTRIBUTIONS_MAX];
extern int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];


// declarations of all the functions used

int min (int a, int b);
int max (int a, int b);
double myClock ();

short int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net);
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net);

void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[]);
void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[]);
void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z);
double lbmStress (double f[]);
void lbmCalculateBC (double f[], unsigned int site_data,double  *vx, double *vy, double *vz, LBM *lbm);
void lbmInit (char *system_file_name, char *parameters_file_name, char *checkpoint_file_name,
	      LBM *lbm, Net *net);
int lbmCycle (int write_checkpoint, int check_convergence, int *is_converged, LBM *lbm, Net *net);
void lbmEnd (LBM *lbm);

void netFindTopology (Net *net);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

void lbmReadAndSetConfig (LBM *lbm, Net *net);
void lbmReadParameters (LBM *lbm);
void lbmSetInitialConditions (LBM *lbm, Net *net);
void lbmWriteConfig (int stability, char *output_file_name, int is_checkpoint, LBM *lbm, Net *net);


#endif                  // __config_h__
