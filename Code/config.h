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


#include <pthread.h>


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
#define COMMS_LEVELS                   2
#define COLLISION_TYPES                6
#define MACHINES_MAX                   4

#define PIXELS_X                       1024
#define PIXELS_Y                       1024
#define COLOURED_PIXELS_PER_PROC_MAX   PIXELS_X * PIXELS_Y
#define IMAGE_SIZE                     PIXELS_X * PIXELS_Y
#define COLOURED_PIXELS_PER_PROC_MAX   PIXELS_X * PIXELS_Y

#define VIS_SIMD_SIZE                  3


#define PixelI(i)   i >> 16
#define PixelJ(i)   i & 65535


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


extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_cond_t network_send_frame;

extern int send_array_length;


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
  int    period;
  
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


struct ProcBlock
{
  int *proc_id;
};


struct BlockLocation
{
  short int i, j, k;
};


struct SiteLocation
{
  short int i, j, k;
};


struct LBM
{
  char *system_file_name;
  
  double tau;
  double viscosity;
  double omega;
  double lattice_to_system;
  double tolerance;
  double conv_error;
  
  int total_fluid_sites;
  int site_min_x, site_min_y, site_min_z;
  int site_max_x, site_max_y, site_max_z;
  int inlets, outlets;
  int cycles_max;
  int period;
  int conv_freq;
  
  float *block_density;
  
  int *block_map;
  
  short int *fluid_sites_per_block;
};


struct NeighProc
{
  int id;
  int fs;
  
  short int *f_data;
  
  int f_head;
  int *f_recv_iv;
};


struct Net
{
  int id;
  int procs;
  int neigh_procs;
  int err;
  int my_inter_sites, my_inner_sites;
  int my_inner_collisions[COLLISION_TYPES];
  int my_inter_collisions[COLLISION_TYPES];
  int my_sites;
  int shared_fs;
  
  int *machine_id;
  int *procs_per_machine;
  int *fluid_sites;
  
  short int *from_proc_id_to_neigh_proc_index;
  short int *proc_id;
  short int *cluster_id;
  
  DataBlock *data_block;
  DataBlock *map_block;
  
  ProcBlock *proc_block;
  
  NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
#ifndef NOMPI
  MPI_Status status[4];
  
  MPI_Request **req;
#endif
  double dd_time, bm_time, fr_time, fo_time;
};


struct ColPixel
{
  union {float r, v;};
  union {float g, t;};
  float b;
  
  int i;
};


struct Cluster
{
  float minmax_x[2], minmax_y[2], minmax_z[2];
  
  float x[4];
  
  unsigned short int blocks_x, blocks_y, blocks_z;
  unsigned short int block_min[3];
};


struct Vis
{
  char *image_file_name;
  
  float half_dim[VIS_SIMD_SIZE];
  float system_size;
};


struct Screen
{
  float vtx[3];
  float dir1[3];
  float dir2[3];
  
  float max_x, max_y;
  float scale_x, scale_y;
  
  float zoom;
  
  int pixels_x, pixels_y;
};


struct Viewpoint
{
  float x[VIS_SIMD_SIZE];
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

extern Cluster *cluster;

extern float **cluster_voxel;

extern float ***cluster_flow_field;


extern short int *f_data;

extern int *f_recv_iv;

extern unsigned int *net_site_data;

extern double *inlet_density;
extern double *inlet_density_avg, *inlet_density_amp, *inlet_density_phs;
extern double *outlet_density;
extern double *outlet_density_avg, *outlet_density_amp, *outlet_density_phs;


extern int col_pixels, col_pixels_max;
extern int col_pixels_recv[ MACHINES_MAX-1 ];
//extern int col_pixels_locked;

extern int *col_pixel_id;

// extern ColPixel *col_pixel_send;
extern ColPixel col_pixel_send[ (MACHINES_MAX-1)*COLOURED_PIXELS_PER_PROC_MAX ];
extern ColPixel *col_pixel_recv;

//extern ColPixel *col_pixel_locked;


extern int is_bench;


extern int sites_x, sites_y, sites_z;
extern int blocks_x, blocks_y, blocks_z;
extern int blocks_yz, blocks;
extern int block_size, block_size2, block_size3, block_size_1;
extern int shift;
extern int sites_in_a_block;

extern double lbm_stress_par;
extern double lbm_density_min, lbm_density_max;
extern double lbm_velocity_min, lbm_velocity_max;
extern double lbm_stress_min, lbm_stress_max;

extern int net_machines;

extern int vis_mode;
extern int vis_flow_field_type;
extern int vis_image_freq;
extern int vis_pixels_max;
extern int vis_compositing;

extern float block_size_f;
extern float block_size_inv;
extern float vis_flow_field_value_max_inv;
extern float vis_absorption_factor;
extern float vis_cutoff;
extern float vis_flow_field_cutoff;
extern float vis_t_min;
extern float vis_value;


extern int cluster_blocks_vec[3];
extern int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;

extern float block_size_f;
extern float block_size_inv;
extern float vis_flow_field_value_max_inv;
extern float vis_cutoff;
extern float vis_flow_field_cutoff;
extern float vis_t_min;


extern float ray_dir[3];
extern float ray_inv[3];
extern float ray_col[3];

extern int clusters;


extern Screen screen;

extern Viewpoint viewpoint;

extern Vis vis;


// declarations of all the functions used

int min (int a, int b);
int max (int a, int b);
int nint (float a);
double myClock ();

int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net);
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net);

void lbmFeq (double f[], double *density, double *v_x, double *v_y, double *v_z, double f_eq[]);
void lbmFeq (double density, double v_x, double v_y, double v_z, double f_eq[]);
void lbmDensityAndVelocity (double f[], double *density, double *v_x, double *v_y, double *v_z);
void lbmStress (double f[], double *stress);
void lbmCalculateBC (double f[], unsigned int site_data, double *density, double *vx, double *vy, double *vz, double f_neq[]);
int lbmCollisionType (unsigned int site_data);
void lbmInit (char *system_file_name, LBM *lbm, Net *net);
void lbmSetInitialConditions (Net *net);
void lbmUpdateFlowField (int perform_rt, int i, double density, double vx, double vy, double vz, double f_neq[]);
int lbmCycle (int cycle_id, int time_step, int perform_rt, LBM *lbm, Net *net);
void lbmEnd (LBM *lbm);

int netFindTopology (Net *net, int *depths);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

void lbmReadConfig (LBM *lbm, Net *net);

#ifdef STEER
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net, SteerParams *steer);
void lbmUpdateParameters (LBM *lbm, SteerParams *steer);
#else
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net);
#endif

void lbmWriteConfig (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmVaryBoundaryDensities (int cycle_id, int time_step, LBM *lbm);

void rtInit (Net *net);
void rtUpdateColour (float dt, float palette[], float col[]);
void rtRayTracingVR (void (*ColourPalette) (float value, float col[]));
void rtRayTracingIS (void (*ColourPalette) (float value, float col[]));
void rtEnd (void);


void visProject (float p1[], float p2[]);
void visWritePixel (ColPixel *col_pixel);
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
void visInit (Net *net, Vis *vis);
void visRenderA (void (*ColourPalette) (float value, float col[]), Net *net);
void visRenderB (int stream_image, int write_image, char *image_file_name, void (*ColourPalette) (float value, float col[]),
		 Net *net);

#ifdef STEER
void visReadParameters (char *parameters_file_name, Net *net, Vis *vis, SteerParams *steer);
void visUpdateParameters (LBM *lbm, Vis *vis, SteerParams *steer);
#else
void visReadParameters (char *parameters_file_name, Net *net, Vis *vis);
#endif

void visEnd (void);

#endif                  // __config_h__

