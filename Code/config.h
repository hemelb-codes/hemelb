// In this file all the declarations and some simple functions are
// reported.
// The structs are defined here.

// Global coordinate means coordinate within the entire system, not the
// coordinate on one proc.
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

#include "colpixel.h"
#include "configconstants.h"
#include "lb.h"

#ifndef NO_STEER
#include <pthread.h>
#include <semaphore.h>
#endif


#define PI           3.14159265358979323846264338327950288
#define DEG_TO_RAD   PI / 180.0
#define EPSILON      1.0e-30

#define VON_MISES_STRESS   +1.0
#define SHEAR_STRESS       -1.0

#define STABLE                 1
#define UNSTABLE               0
#define STABLE_AND_CONVERGED   2

#define MACROSCOPIC_PARS   5
#define DENSITY            0
#define VELOCITY           1
#define STRESS             2

#define COMMS_LEVELS                   2


#define RECV_BUFFER_A   0
#define RECV_BUFFER_B   1

#define STEERABLE_PARAMETERS   20

#define REFERENCE_PRESSURE             80.0           // 80 mmHg
#define mmHg_TO_PASCAL                 133.3223874
#define BLOOD_DENSITY                  1000.0        // 1000 Kg m^(-3)
#define BLOOD_VISCOSITY                0.004         // 0.004 Pascal s
#define PULSATILE_PERIOD               0.857142857   // period of oscillation (in s) is
					             // chosen to be 1 min / 70
					             // beats per min
#define TOL                            1.0e-6

#define VIS_FIELDS                     3

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

// Some sort of coordinates.
struct BlockLocation
{
  short int i, j, k;
};



// Some sort of coordinates.
struct SiteLocation
{
  short int i, j, k;
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

extern double *net_site_nor;
extern unsigned int *net_site_data;

extern double *inlet_density;
extern double *inlet_density_avg, *inlet_density_amp, *inlet_density_phs;
extern double *outlet_density;
extern double *outlet_density_avg, *outlet_density_amp, *outlet_density_phs;


extern int col_pixels, col_pixels_max;
extern int col_pixels_recv[2];

extern int *col_pixel_id;
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
int *netProcIdPointer (int site_i, int site_j, int site_k, Net *net);
unsigned int *netSiteMapPointer (int site_i, int site_j, int site_k, Net *net);



int netFindTopology (Net *net, int *depths);
void netInit (LBM *lbm, Net *net);
void netEnd (Net *net);

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
