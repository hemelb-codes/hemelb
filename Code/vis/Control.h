#ifndef HEMELB_VIS_CONTROL_H
#define HEMELB_VIS_CONTROL_H

#include <vector>

#include "net.h"
#include "lb.h"

#include "vis/ColPixel.h"
#include "vis/GlyphDrawer.h"
#include "vis/StreaklineDrawer.h"
namespace hemelb {
  namespace vis {
    /* the last three digits of the pixel identifier are used to
       indicate if the pixel is coloured via the ray tracing technique
       and/or a glyph and/or a particle/pathlet */
    

    // Class to control and use the effects of different visualisation
    // methods.
    class Control
    {
    public:
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


      struct Vis
      {
	// As far as I can tell image_file_name member is never used
	char *image_file_name;
  
	float half_dim[3];
	float system_size;
      };

      Control();
      ~Control();
      void initLayers(Net* net);
      /* This next bit was Hywel's incomplete attempt to refactor the
	 vis Control subsystem. I'm going to remove it for now and
	 reintegrate it later. */

      /* // Adds a layer to the visualisation. Note that rendering is done */
      /* // in the order in which layers are added. */
      /* void addLayer(Layer *newLayer); */
    
      /* // Combines the output of all visualisation methods in the order */
      /* // added. */
      /* void render(); */

      
      
      void project (float p1[], float p2[]);
      void writePixel (ColPixel *col_pixel);
      void mergePixels (ColPixel *col_pixel1, ColPixel *col_pixel2);
      void rotate (float sin_1, float cos_1,
		   float sin_2, float cos_2,
		   float  x1, float  y1, float  z1,
		   float *x2, float *y2, float *z2);
      
      void setProjection (int pixels_x, int pixels_y,
			  float longitude, float latitude,
			  float zoom);
      
      void renderLine (float x1[], float x2[]);
      

      void streaklines(int time_step, int period, Net *net);
      void restart();

      void updateImageSize (int pixels_x, int pixels_y);
      void render (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net);
      void writeImage (int recv_buffer_id, char *image_file_name,
		       void (*ColourPalette) (float value, float col[]));
      void readParameters (char *parameters_file_name, LBM *lbm, Net *net);
      void calculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm);
      void compositeImage (int recv_buffer_id, Net *net);

      // better public member vars than globals!
      int mode;
      int image_freq;
      int shouldDrawStreaklines;

      float velocity_threshold_max_inv;
      float stress_threshold_max_inv;
      float density_threshold_min, density_threshold_minmax_inv;
      float ctr_x, ctr_y, ctr_z;
      float streaklines_per_pulsatile_period, streakline_length;
      double mouse_pressure, mouse_stress;
      float brightness;

      Screen screen;
      Viewpoint viewpoint;
      
      int col_pixels_recv[2]; // number received?
      //ColPixel col_pixel_recv[2][COLOURED_PIXELS_MAX];
      ColPixel* col_pixel_recv[2];

      float physical_pressure_threshold_min;
      float physical_pressure_threshold_max;
      float physical_velocity_threshold_max;
      float physical_stress_threshold_max;
      
      int mouse_x, mouse_y;
      
    private:
      /* Hywel's partial refactor */
      /* int nLayers; */
      /* std::vector<Layer *> myLayers; */
      
      int col_pixels, // number of ColPixels (?)
	col_pixels_max; // max permitted of above
      
      int *col_pixel_id; // array of pixel IDs
      
      // These don't belong in vis:: namespace
      /*double pressure_min, pressure_max;
	double velocity_min, velocity_max;
	double stress_min, stress_max;*/
      //double time;

      int time_step, cycle;
      int period, inlets;
      int pixels_max;


      int perform_rendering;

      ColPixel col_pixel_send[COLOURED_PIXELS_MAX];

      Vis vis;

      GlyphDrawer *myGlypher;
      StreaklineDrawer *myStreaker;

    };
    
    extern Control *controller;
  }
}

#endif // HEMELB_VIS_CONTROL_H
