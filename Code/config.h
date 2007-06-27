// In this file all the declarations and some simple functions are reported.

#ifndef __config_h__
#define __config_h__

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <rpc/xdr.h>

int STABLE   = 1;
int UNSTABLE = 0;

// the constants needed to define the configuration of the lattice
// sites follow

unsigned int SOLID_TYPE  = 0U;
unsigned int FLUID_TYPE  = 1U;
unsigned int INLET_TYPE  = 2U;
unsigned int OUTLET_TYPE = 3U;
unsigned int NULL_TYPE   = 4U;

unsigned int BOUNDARIES              = 4U;
unsigned int INLET_BOUNDARY          = 0U;
unsigned int OUTLET_BOUNDARY         = 1U;
unsigned int WALL_BOUNDARY           = 2U;
unsigned int CHARACTERISTIC_BOUNDARY = 3U;

unsigned int SITE_TYPE_BITS       = 2U;
unsigned int BOUNDARY_CONFIG_BITS = 14U;
unsigned int BOUNDARY_DIR_BITS    = 4U;
unsigned int BOUNDARY_ID_BITS     = 10U;

unsigned int BOUNDARY_CONFIG_SHIFT = 2U;   // SITE_TYPE_BITS;
unsigned int BOUNDARY_DIR_SHIFT    = 16U;  // BOUNDARY_CONFIG_SHIFT + BOUNDARY_CONFIG_BITS;
unsigned int BOUNDARY_ID_SHIFT     = 20U;  // BOUNDARY_DIR_SHIFT + BOUNDARY_DIR_BITS;

unsigned int SITE_TYPE_MASK       = ((1U <<  2U) - 1U);         // ((1U << SITE_TYPE_BITS) - 1U);
unsigned int BOUNDARY_CONFIG_MASK = ((1U << 14U) - 1U) << 2U;   // ((1U << BOUNDARY_CONFIG_BITS) - 1U) << BOUNDARY_CONFIG_SHIFT;
unsigned int BOUNDARY_DIR_MASK    = ((1U <<  4U) - 1U) << 16U;  //((1U << BOUNDARY_DIR_BITS) - 1U)    << BOUNDARY_DIR_SHIFT;
unsigned int BOUNDARY_ID_MASK     = ((1U << 10U) - 1U) << 20U;  // ((1U << BOUNDARY_ID_BITS) - 1U)     << BOUNDARY_ID_SHIFT
unsigned int CHARACTERISTIC_MASK  = 1U << 31U;

// square of the speed of sound

double Cs2 = 3.0 / 8.0;


// parameters related to the lattice directions

int e_x[] = { 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};
int e_y[] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1,-1, 1,-1, 1};
int e_z[] = { 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1};
int inv_dir[] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};


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


double *f_old, *f_new;

int *f_id;

Velocity *vel;


short int f_data[4*SHARED_DISTRIBUTIONS_MAX];

double f_to_send[SHARED_DISTRIBUTIONS_MAX];
double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

int f_send_id[SHARED_DISTRIBUTIONS_MAX];
int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];


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
void CalculateBC (double f[], unsigned int site_data,double  *vx, double *vy, double *vz, LBM *lbm);

void netFindTopology (Net *net);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

void lbmReadAndSetConfig (LBM *lbm, Net *net);
void lbmReadParameters (LBM *lbm);
void lbmSetInitialConditions (LBM *lbm, Net *net);
void lbmWriteConfig (int stability, char *output_file_name, int is_checkpoint, LBM *lbm, Net *net);



// some simple functions

int min (int a, int b)
{
  if (a < b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


int max (int a, int b)
{
  if (a > b)
    {
      return a;
    }
  else
    {
      return b;
    }
}


double myClock ()
{
  return (double)clock () / (double)CLOCKS_PER_SEC;
}


#endif                  // __config_h__
