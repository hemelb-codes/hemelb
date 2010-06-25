#ifndef __rt_h_
#define __rt_h_

#include "net.h"
#include "lb.h"

#define RT               (1 << 28)
#define GLYPH            (1 << 29)
#define STREAKLINE       (1 << 30)
#define PIXEL_ID_MASK    ((1 << 28) - 1)

#define PixelI(i)        ((i >> 14) & 8191)
#define PixelJ(i)        (i & 16383)
#define PixelId(i,j)     ((i << 14) | j)

extern float **cluster_voxel;

extern float ***cluster_flow_field;

extern int col_pixels, col_pixels_max;
extern int col_pixels_recv[2];

extern int *col_pixel_id;

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

extern int clusters;


extern int glyphs;






struct AABB
{
  float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
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

extern ColPixel col_pixel_send[COLOURED_PIXELS_MAX];
extern ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];


// Some sort of coordinates.
struct BlockLocation
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

void rawWritePixel(ColPixel*, unsigned int*, unsigned char [],
		   void (*ColourPalette) (float, float []));
void makePixelColour(unsigned char& red, unsigned char& green, unsigned char& blue,
  int rawRed, int rawGreen, int rawBlue);

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


struct Vis
{
  char *image_file_name;
  
  float half_dim[3];
  float system_size;
};


extern Vis vis;
extern Glyph *glyph;

void visProject (float p1[], float p2[]);
void visWritePixel (ColPixel *col_pixel);
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



extern float ray_dir[3];
extern float ray_inv[3];
extern float ray_vel_col[3];
extern float ray_stress_col[3];
extern float ray_length;
extern float ray_t_min;
extern float ray_density;
extern float ray_stress;
extern Cluster *cluster;
extern Screen screen;
extern Viewpoint viewpoint;


#endif // __rt_h_
