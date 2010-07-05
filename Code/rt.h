#ifndef __rt_h_
#define __rt_h_

#include "net.h"
#include "lb.h"

#include "glyphDrawer.h"
#include "colpixel.h"
#include "rayTracer.h"

#ifndef NO_STREAKLINES
  #include "streaklineDrawer.h"
#endif // NO_STREAKLINES


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
extern float vis_streaklines_per_pulsatile_period, vis_streakline_length;

extern int vis_mouse_x, vis_mouse_y;
extern int vis_perform_rendering;
extern int vis_mode;

extern ColPixel col_pixel_send[COLOURED_PIXELS_MAX];
extern ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];


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



void rawWritePixel(ColPixel*, unsigned int*, unsigned char [],
		   void (*ColourPalette) (float, float []));
void makePixelColour(unsigned char& red, unsigned char& green, unsigned char& blue,
  int rawRed, int rawGreen, int rawBlue);

struct Vis
{
  char *image_file_name;
  
  float half_dim[3];
  float system_size;
};

extern Vis vis;

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
void visInit (Net *net, Vis *vis);

void visStreaklines(int time_step, int period, Net *net);
void visRestart();

void visUpdateImageSize (int pixels_x, int pixels_y);
void visRender (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net);
void visWriteImage (int recv_buffer_id, char *image_file_name,
		    void (*ColourPalette) (float value, float col[]));
void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis);
void visCalculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm);
void visEnd ();

extern Screen screen;
extern Viewpoint viewpoint;


#endif // __rt_h_
