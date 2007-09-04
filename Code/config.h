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

#ifdef RG
#include <pthread.h>
#endif // RG

#define MACROSCOPIC_PARS   5
#define DENSITY            0
#define VELOCITY           1
#define STRESS             2

#define NEIGHBOUR_PROCS_MAX            52
#define SHARED_DISTRIBUTIONS_MAX       4000000
#define COLOURED_PIXELS_PER_PROC_MAX   512 * 512


extern float EPSILON;

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

extern int e_x[15];
extern int e_y[15];
extern int e_z[15];
extern int inv_dir[15];


// data structures useful to define the simulation set-up, construct
// the system and partition it

#ifdef RG

extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_cond_t network_send_frame;

struct Pixel {
	int x, y;
	float r, g, b;
};

// extern Pixel* send_array;
extern int send_array_length;
extern int send_frame_count;

extern int compressedFrameSize;

extern unsigned char *pixelData;
extern unsigned char *compressedData;

#endif // RG


struct DataBlock
{
  unsigned int *site_data;
};


struct BlockLocation
{
  short int i, j, k;
};


struct SiteLocation ///
{
  short int i, j, k;
};


struct LBM
{
  char *system_file_name;
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
  int inlet_sites; ///
  
  float *block_density;
  
  int *block_map;
  
  short int *fluid_sites_per_block;
  
  SiteLocation *site_location_a; ///
  SiteLocation *site_location_b; ///
  
  DataBlock *iter_block; ///
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
  
  double **d_to_send_p;
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
  int convergence_count;
  
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
  
  double timing[8];
};


struct BlockData
{
  DataBlock *p;
  
  float t;
  float min_x, min_y, min_z;
  
  int i, j, k;
};


struct Cluster
{
  float block_min_x, block_min_y, block_min_z;
  
  int blocks_yz, blocks;
  
  unsigned short int blocks_x, blocks_y, blocks_z;
  unsigned short int block_min_i, block_min_j, block_min_k;
};


struct RT
{
  char *image_file_name;
  
  int image_frequency;
  int flow_field_type;
  int is_isosurface;
  int blocks_x, blocks_y, blocks_z;
  int blocks_yz, blocks;
  int block_size, block_size2, block_size3, block_size_1;
  int shift, shift2;
  int coloured_pixels_max;
  int coloured_pixels;
  int *coloured_pixel_id;
  int coloured_pixel_ids_max;
  int ray_tracing_count;
  
  short int clusters;
  
  float flow_field_value_max_inv[3];
  float absorption_factor;
  float cutoff;
  float t_min;
  float block_size_inv;
  float system_size;
  float half_dim_x, half_dim_y, half_dim_z;
  
  float coloured_pixel_send[6*COLOURED_PIXELS_PER_PROC_MAX];
  float *coloured_pixel_recv;
  
  Cluster *cluster;
};


struct Screen
{
  float ctr_x, ctr_y, ctr_z;
  float max_x, max_y;
  
  float dir1x, dir1y, dir1z;
  float dir2x, dir2y, dir2z;
  float par_x, par_y, par_z;
  
  float zoom, dist;
  
  int pixels_x, pixels_y;
};


struct Viewpoint
{
  float pos_x, pos_y, pos_z;
  float sin_1, cos_1;
  float sin_2, cos_2;
};


struct Ray
{
  float dir_x, dir_y, dir_z;
  float inv_x, inv_y, inv_z;
  float col_r, col_g, col_b;
};


struct AABB
{
  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
};


extern double *f_old, *f_new;

extern int *f_id;

extern Velocity *vel;

extern float *flow_field;

extern double *d;

extern double **nd_p;


extern short int f_data[4*SHARED_DISTRIBUTIONS_MAX];

extern double f_to_send[SHARED_DISTRIBUTIONS_MAX];
extern double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

extern int f_send_id[SHARED_DISTRIBUTIONS_MAX];
extern int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];


extern Screen screen;

extern Viewpoint viewpoint;

extern Ray ray;


// declarations of all the functions used

int min (int a, int b);
int max (int a, int b);
double myClock ();

short int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net);
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net);

void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[]);
void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[]);
void lbmVelocity (double f[], double *v_x, double *v_y, double *v_z);
void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z);
double lbmStress (double f[]);
void lbmCalculateBC (double f[], unsigned int site_data, double *density, double *vx, double *vy, double *vz, LBM *lbm);
void lbmInit (char *system_file_name, char *checkpoint_file_name, LBM *lbm, Net *net);
void lbmSetOptimizedInitialConditions (LBM *lbm, Net *net);
int lbmCycle (int write_checkpoint, int check_convergence, int perform_rt, int *is_converged, LBM *lbm, Net *net);
void lbmEnd (LBM *lbm);

void netFindTopology (Net *net);
void netInit (LBM *lbm, Net *net, RT *rt);
void netEnd (Net *net, RT *rt);

void lbmReadConfig (LBM *lbm, Net *net);
void lbmReadParameters (char * parameters_file_name, LBM *lbm, Net *net);
void lbmSetInitialConditions (LBM *lbm, Net *net);
void lbmWriteConfig (int stability, char *output_file_name, int is_checkpoint, LBM *lbm, Net *net);


void rtProject (float px1, float py1, float pz1, float *px2, float *py2);
void rtRayTracing (void (*AbsorptionCoefficients) (float flow_field_data, float t1, float t2,
						   float cutoff, float *r, float *g, float *b),
		   LBM *lbm, Net *net, RT *rt);
void rtRotate (float sin_1, float cos_1,
	       float sin_2, float cos_2,
	       float  x1, float  y1, float  z1,
	       float *x2, float *y2, float *z2);
void rtProjection (float ortho_x, float ortho_y,
		   int pixels_x, int pixels_y,
		   float ctr_x, float ctr_y, float ctr_z,
		   float rad,
		   float longitude, float latitude,
		   float dist,
		   float zoom);
void rtReadParameters (char * parameters_file_name, RT *rt, Net * net);
void rtInit (char *image_file_name, RT *rt);
void rtEnd (RT *rt);

#endif                  // __config_h__
