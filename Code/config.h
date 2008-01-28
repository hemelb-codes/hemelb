// In this file all the declarations and some simple functions are reported.

#ifndef __config_h__
#define __config_h__

#ifndef NOMPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef int64_t
#define int64_t long int
#endif

#ifdef XT3
#include "types.h"
#include "xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif


#ifdef BENCH
#undef RG
#undef STEER
#endif


#ifdef RG
#include <pthread.h>
#endif // RG


#ifdef STEER
#include "ReG_Steer_types.h"
#include "ReG_Steer_Appside.h"
#endif // STEER


#define EPSILON    1.e-30

#define STABLE     1
#define UNSTABLE   0

#define MACROSCOPIC_PARS   5
#define DENSITY            0
#define VELOCITY           1
#define STRESS             2

#define NEIGHBOUR_PROCS_MAX            52
#define SHARED_DISTRIBUTIONS_MAX       4000000
#define COMMS_LEVELS                   2
#define COLLISION_TYPES                6
#define MACHINES_MAX                   4
#define PROCS_MAX                      8 * 1024

#define STREAMLINES_MAX                1000

#define PIXELS_X                       1024
#define PIXELS_Y                       1024
#define COLOURED_PIXELS_PER_PROC_MAX   PIXELS_X * PIXELS_Y
#define IMAGE_SIZE                     PIXELS_X * PIXELS_Y
#define SIMD_SIZE                      2
#define VIS_VEC_SIZE                   4


extern double PI;

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

#ifndef NOMPI
extern MPI_Datatype MPI_col_pixel_type;
#endif

// data structures useful to define the simulation set-up, construct
// the system and partition it

#ifdef RG

extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_cond_t network_send_frame;

extern int send_array_length;

#endif // RG


#ifdef STEER

// this is here so that I can transfer all params and data in one
// chunk in one MPI_Bcast rather than loads of separate ones...
struct SteerParams
{
  // reg
  int    status;
  int    num_recvd_cmds;
  int    recvd_cmds[REG_MAX_NUM_STR_CMDS];
  int    num_params_changed;
  
  // lbm
  double tau;
  double tolerance;
  int    max_cycles;
  int    conv_freq;
#ifndef TD
  int    check_freq;
#else
  int    period;
#endif
  
  // vis
  float  ctr_x, ctr_y, ctr_z;
  float  longitude;
  float  latitude;
  float  zoom;
  int    image_freq;
  int    flow_field_type;
  int    mode;
  float  abs_factor;
  float  cutoff;
  float  max_density;
  float  max_velocity;
  float  max_stress;
};

#endif // STEER


struct DataBlock
{
  unsigned int *site_data;
};


struct BlockLocation
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
  double conv_error;
  double density_min, density_max;
  double velocity_min, velocity_max;
  double stress_min, stress_max;
  
  int flow_field_type;
  int total_fluid_sites;
  int site_min_x, site_min_y, site_min_z;
  int site_max_x, site_max_y, site_max_z;
  int inlets, outlets;
  int cycles_max;
#ifndef TD
  int checkpoint_freq;
#endif
  int period;
  int conv_freq;
  int is_checkpoint;
  
  float *block_density;
  
  int *block_map;
  
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
  
#ifndef BENCH
  double **d_to_send_p;
#endif
};


struct Net
{
  int id;
  int procs, machines;
  int neigh_procs;
  int err;
  int my_inter_sites, my_inner_sites;
  int my_inner_collisions[COLLISION_TYPES];
  int my_inter_collisions[COLLISION_TYPES];
  int my_inner_collisions_simd[COLLISION_TYPES];
  int my_inter_collisions_simd[COLLISION_TYPES];
  int my_sites;
  int shared_fs;
  
  int *machine_id;
  int *procs_per_machine;
  
  short int *from_proc_id_to_neigh_proc_index;
  short int *proc_id;
  
  unsigned int *site_data;
  
  DataBlock *data_block;
  DataBlock *map_block;
  
  NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
#ifndef NOMPI
  MPI_Status status[4];
  
  MPI_Request **req;
#endif
  double dd_time, bm_time, fr_time, fo_time;
};


struct ColPixel
{
  float r, g, b;
  float t;
  
  short int i, j;
};


struct Cluster
{
  float x[VIS_VEC_SIZE];
  
  unsigned short int blocks_x, blocks_y, blocks_z;
  unsigned short int block_min[3];
};


struct Vis
{
  float flow_field_value_max_inv;
  float cutoff;
  float t_min;
  
  int mode;
  
  char *image_file_name;
  
  int image_freq;
  int col_pixels, col_pixels_max, col_pixels_locked;
  int col_pixels_recv[ MACHINES_MAX-1 ];
  int *col_pixel_id;
  int pixels_max;
  int seeds, seeds_max;
  int streamlines, streamlines_max;
  
  float velocity_max_inv;
  float absorption_factor;
  float half_dim[VIS_VEC_SIZE];
  float system_size;
  
  double stress_par;
  
  float *seed;
  float *streamline;
  
  ColPixel col_pixel_send[ (MACHINES_MAX-1)*COLOURED_PIXELS_PER_PROC_MAX ];
  ColPixel *col_pixel_recv;
  ColPixel *col_pixel_locked;
};


struct Screen
{
  float vtx[VIS_VEC_SIZE];
  float dir1[VIS_VEC_SIZE];
  float dir2[VIS_VEC_SIZE];
  
  float max_x, max_y;
  float scale_x, scale_y;
  
  float zoom;
  
  int pixels_x, pixels_y;
};


struct Viewpoint
{
  float x[VIS_VEC_SIZE];
  float sin_1, cos_1;
  float sin_2, cos_2;
  float dist;
};


struct AABB
{
  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
};


extern double *f_old, *f_new;

extern int *f_id;

#ifndef TD
extern double *vel;
#endif

extern float *flow_field;

#ifndef BENCH
extern double *d;

extern double **nd_p;
#endif

extern Cluster *cluster;


extern short int f_data[4*SHARED_DISTRIBUTIONS_MAX];
extern double f_to_send[SHARED_DISTRIBUTIONS_MAX];
extern double f_to_recv[SHARED_DISTRIBUTIONS_MAX];

extern int f_send_id[SHARED_DISTRIBUTIONS_MAX];
extern int f_recv_iv[SHARED_DISTRIBUTIONS_MAX];

extern float streamline_to_send[NEIGHBOUR_PROCS_MAX][VIS_VEC_SIZE*STREAMLINES_MAX];
extern float streamline_to_recv[NEIGHBOUR_PROCS_MAX][1+VIS_VEC_SIZE*STREAMLINES_MAX];

extern int streamlines_to_send[NEIGHBOUR_PROCS_MAX];
extern int streamlines_to_recv[NEIGHBOUR_PROCS_MAX];


extern int sites_x, sites_y, sites_z;
extern int blocks_x, blocks_y, blocks_z;
extern int blocks_yz, blocks;
extern int block_size, block_size2, block_size3, block_size_1;
extern int shift;
extern int sites_in_a_block;

extern float block_size_inv;


extern float ray_dir[VIS_VEC_SIZE];
extern float ray_inv[VIS_VEC_SIZE];
extern float ray_col[VIS_VEC_SIZE];

extern short int clusters;


extern Screen screen;

extern Viewpoint viewpoint;

extern Vis vis;


// declarations of all the functions used

int min (int a, int b);
int max (int a, int b);
int nint (float a);
double myClock ();

short int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net);
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net);

void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[]);
void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[]);
void lbmVelocity (double f[], double *v_x, double *v_y, double *v_z);
void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z);
void lbmStressSIMD (double f[], double stress[], LBM *lbm);
void lbmStress (double f[], double *stress, LBM *lbm);
void lbmCalculateBC (double f[], unsigned int site_data, double *density, double *vx, double *vy, double *vz, LBM *lbm);
int lbmCollisionType (unsigned int site_data);
void lbmInit (char *system_file_name, char *checkpoint_file_name, LBM *lbm, Net *net);
void lbmSetInitialConditions (LBM *lbm, Net *net);
void lbmUpdateFlowFieldSIMD (int i, double f_neq[], double density[], double vx[], double vy[], double vz[], LBM *lbm);
void lbmUpdateFlowField (int i, double f_neq[], double density, double vx, double vy, double vz, LBM *lbm);
#ifndef TD
int lbmCycle (int write_checkpoint, int check_conv, int perform_rt, int *is_converged, LBM *lbm, Net *net);
#else
int lbmCycle (int cycle_id, int time_step, int check_conv, int perform_rt, int *is_converged, LBM *lbm, Net *net);
#endif
void lbmEnd (LBM *lbm);

int netFindTopology (Net *net, int *depths);
void netInit (LBM *lbm, Net *net, int proc_sites[]);
void netEnd (Net *net);

void lbmReadConfig (LBM *lbm, Net *net);

#ifdef STEER
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net, SteerParams *steer);
#else
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net);
#endif

#ifdef STEER
void lbmUpdateParameters (LBM *lbm, SteerParams *steer);
#endif

void lbmWriteConfig (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmSetInitialConditionsWithCheckpoint (LBM *lbm, Net *net);

void rtInit (Net *net, Vis *vis);
void rtUpdateColour (float dt, float palette[], float col[]);
void rtRayTracing (void (*ColourPalette) (float value, float col[]), Net *net, Vis *vis);
void rtEnd (Vis *vis);

void slInit (Net *net, Vis *vis);
void slStreamlines (Net *net, Vis *vis);
void slEnd (Vis *vis);


void visProject (float p1[], float p2[]);
void visWritePixel (float col[], float t, int i, int j, Vis *vis);
void visRotate (float sin_1, float cos_1,
		float sin_2, float cos_2,
		float  x1, float  y1, float  z1,
		float *x2, float *y2, float *z2);
void visProjection (float ortho_x, float ortho_y,
		    int pixels_x, int pixels_y,
		    float ctr_x, float ctr_y, float ctr_z,
		    float rad,
		    float longitude, float latitude,
		    float dist,
		    float zoom);
void visInit (char *image_file_name, Net *net, Vis *vis);
void visRenderA (void (*ColourPalette) (float value, float col[]), Net *net, Vis *vis);
void visRenderB (Net *net, Vis *vis);

#ifdef STEER
void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis, SteerParams *steer);
#else
void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis);
#endif

#ifdef STEER
void visUpdateParameters (LBM *lbm, Vis *vis, SteerParams *steer);
#endif

void visEnd (Net *net, Vis *vis);

#endif                  // __config_h__

