// In this file all the declarations and some simple functions are reported.

#ifndef __config_h__
#define __config_h__

#ifndef NOMPI
#ifdef XT3
#include <mpi.h>
#else
#include "mpi.h"
#endif
#endif

#ifndef int64_t
#define int64_t long int
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#ifndef NO_STEER
#include <pthread.h>
#include <semaphore.h>
#endif


#define PI           3.14159265358979323846264338327950288
#define DEG_TO_RAD   PI / 180.0
#define EPSILON      1.0e-30

#define STABLE                 1
#define UNSTABLE               0
#define STABLE_AND_CONVERGED   2

#define MACROSCOPIC_PARS   5
#define DENSITY            0
#define VELOCITY           1
#define STRESS             2

#define NEIGHBOUR_PROCS_MAX            64
#define COMMS_LEVELS                   2
#define COLLISION_TYPES                6

#define RECV_BUFFER_A   0
#define RECV_BUFFER_B   1

#define COLOURED_PIXELS_MAX    2048 * 2048
#define STEERABLE_PARAMETERS   20

#define REFERENCE_PRESSURE             80.0           // 80 mmHg
#define mmHg_TO_PASCAL                 133.3223874
#define BLOOD_DENSITY                  1000.0        // 1000 Kg m^(-3)
#define BLOOD_VISCOSITY                0.004         // 0.004 Pascal s
#define PULSATILE_PERIOD               0.857142857   // period of oscillation (in s) is
					             // chosen to be 1 min / 70
					             // beats per min
#define TOL                            1.e-6

// the last three digits of the pixel identifier are used to indicate
// if the pixel is coloured via the ray tracing technique and/or a glyph
// and/or a particle/pathlet
#define RT               (1 << 28)
#define GLYPH            (1 << 29)
#define STREAKLINE       (1 << 30)
#define PIXEL_ID_MASK    ((1 << 28) - 1)

#define PixelI(i)        ((i >> 14) & 8191)
#define PixelJ(i)        (i & 16383)
#define PixelId(i,j)     ((i << 14) | j)


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


#ifndef NO_STEER
extern pthread_mutex_t network_buffer_copy_lock;
extern pthread_mutex_t LOCK;
extern pthread_cond_t network_send_frame;

extern sem_t nrl;
extern sem_t connected_sem;
extern sem_t steering_var_lock;

extern bool is_frame_ready;
extern bool connected;
extern bool sending_frame;

extern int send_array_length;
#endif


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
  
  double tau, viscosity;
  double voxel_size;
  double omega;
  double lattice_to_system;
  
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


struct NeighProc
{
  int id;
  int fs;
  
  short int *f_data;
  
  int f_head;
  int *f_recv_iv;
  
  // buffers needed for convergence-enabled simulations
  double *f_to_send;
  double *f_to_recv;
  
  int *f_send_id;
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
  float vel_r, vel_g, vel_b;
  float stress_r, stress_g, stress_b;
  float t, dt;
  float density;
  float stress;
  
  float particle_vel;
  float particle_z;
  
  int particle_inlet_id;
  int i;
};


struct Cluster
{
  float minmax_x[2], minmax_y[2], minmax_z[2];
  
  float x[3];
  
  unsigned short int blocks_x, blocks_y, blocks_z;
  unsigned short int block_min[3];
};


struct Glyph
{
  float x, y, z;
  
  double *f;
};


struct Particle
{
  float x, y, z;
  float vx, vy, vz;
  float vel;
  
  int inlet_id;
};


struct VelSiteData
{
  int proc_id, counter, site_id;
  
  float vx, vy, vz;
};


struct VelocityField
{
  VelSiteData *vel_site_data;
};


struct SL
{
  int counter;
  int particles, particles_max;
  int particle_seeds, particle_seeds_max;
  int particles_to_send_max, particles_to_recv_max;
  int neigh_procs;
  int shared_vs;
  int procs;
  
  VelocityField *velocity_field;
  
  Particle *particle;
  Particle *particle_seed;
  
  float *v_to_send, *v_to_recv;
  
  short int *s_to_send, *s_to_recv;
  
  short int *from_proc_id_to_neigh_proc_index;
  
  struct NeighProc
  {
    int id;
    int send_ps, recv_ps;
    int send_vs, recv_vs;
    
    float *p_to_send, *p_to_recv;
    float *v_to_send, *v_to_recv;
    
    short int *s_to_send, *s_to_recv;
  };
  
  NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];
  
#ifndef NOMPI
  MPI_Status status[4];
  
  MPI_Request *req;
#endif
};


struct Vis
{
  char *image_file_name;
  
  float half_dim[3];
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
  float x[3];
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


// 3 buffers needed for convergence-enabled simulations
extern double *f_to_send;
extern double *f_to_recv;

extern int *f_send_id;

extern int *f_recv_iv;

extern short int *f_data;

extern unsigned int *net_site_data;

extern double *inlet_density;
extern double *inlet_density_avg, *inlet_density_amp, *inlet_density_phs;
extern double *outlet_density;
extern double *outlet_density_avg, *outlet_density_amp, *outlet_density_phs;


extern int col_pixels, col_pixels_max;
extern int col_pixels_recv[2];

extern int *col_pixel_id;

extern ColPixel col_pixel_send[COLOURED_PIXELS_MAX];
extern ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];
extern Glyph *glyph;


extern int is_bench;

// 3 variables needed for convergence-enabled simulations
extern double conv_error;
extern int cycle_tag, check_conv;
extern int is_inlet_normal_available;


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
extern double *lbm_average_inlet_velocity;
extern double *lbm_peak_inlet_velocity;
extern double *lbm_inlet_normal;
extern long int *lbm_inlet_count;

extern int lbm_terminate_simulation;

extern int net_machines;


extern double vis_pressure_min, vis_pressure_max;
extern double vis_velocity_min, vis_velocity_max;
extern double vis_stress_min, vis_stress_max;
extern double vis_time;

extern int vis_time_step, vis_cycle;
extern int vis_period, vis_inlets;
extern int vis_image_freq;
extern int vis_pixels_max;
extern int vis_streaklines;


extern float block_size_f;
extern float block_size_inv;
extern float vis_physical_pressure_threshold_min;
extern float vis_physical_pressure_threshold_max;
extern float vis_physical_velocity_threshold_max;
extern float vis_physical_stress_threshold_max;
extern float vis_density_threshold_min, vis_density_threshold_minmax_inv;
extern float vis_velocity_threshold_max_inv;
extern float vis_stress_threshold_max_inv;
extern float vis_brightness;
extern float vis_ctr_x, vis_ctr_y, vis_ctr_z;
extern double vis_mouse_pressure, vis_mouse_stress;
extern double vis_glyph_length;
extern float vis_streaklines_per_pulsatile_period, vis_streakline_length;

extern int vis_mouse_x, vis_mouse_y;
extern int vis_perform_rendering;
extern int vis_mode;


extern int cluster_blocks_vec[3];
extern int cluster_blocks_z, cluster_blocks_yz, cluster_blocks;


extern float ray_dir[3];
extern float ray_inv[3];
extern float ray_vel_col[3];
extern float ray_stress_col[3];
extern float ray_length;
extern float ray_t_min;
extern float ray_density;
extern float ray_stress;

extern int clusters;


extern int glyphs;


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

int netFindTopology (Net *net, int *depths);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

void lbmReadConfig (LBM *lbm, Net *net);
double lbmCalculateTau (LBM *lbm);
void lbmReadParameters (char *parameters_file_name, LBM *lbm, Net *net);

void lbmWriteConfig (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmWriteConfigASCII (int stability, char *output_file_name, LBM *lbm, Net *net);
void lbmUpdateBoundaryDensities (int cycle_id, int time_step, LBM *lbm);
void lbmUpdateInletVelocities (int time_step, LBM *lbm, Net *net);

void rtInit (Net *net);
void rtUpdateRayData (float *flow_field, float ray_t, float ray_segment, void (*ColourPalette) (float value, float col[]));
void rtRayTracing (void (*ColourPalette) (float value, float col[]));
void rtUpdateClusterVoxel (int i, float density, float velocity, float stress);
void rtEnd (void);


void glyInit (Net *net);
void glyGlyphs (void);
void glyEnd (void);


void slStreakLines (int time_steps, int time_steps_per_cycle, Net *net, SL *sl);
void slRender (int recv_buffer_id, SL *sl);
void slRestart (SL *sl);


void visProject (float p1[], float p2[]);
void visWritePixel (ColPixel *col_pixel);
void rawWritePixel(ColPixel*, unsigned int*, unsigned char [],
		   void (*ColourPalette) (float, float []));
void xdrWritePixel (ColPixel *col_pixel_p, XDR *xdr_p, void (*ColourPalette) (float value, float col[]));
void visMergePixels (ColPixel *col_pixel1, ColPixel *col_pixel2);
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
void visRenderLine (float x1[], float x2[]);
void visInit (Net *net, Vis *vis, SL *sl);
void visUpdateImageSize (int pixels_x, int pixels_y);
void visRender (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net, SL *sl);
void visWriteImage (int recv_buffer_id, char *image_file_name,
		    void (*ColourPalette) (float value, float col[]));
void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis);
void visCalculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm);
void visEnd (SL *sl);

#endif                  // __config_h__
