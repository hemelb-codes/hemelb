#include "config.h"
#include "steering.h"

extern float steer_par[];

void UpdateSteerableParameters (int *vis_perform_rendering, Vis *vis, LBM* lbm)
{
  steer_par[ STEERABLE_PARAMETERS ] = *vis_perform_rendering;

#ifndef NOMPI
    MPI_Bcast (steer_par, STEERABLE_PARAMETERS+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  float longitude, latitude;
  float zoom;
  float velocity_max, stress_max;
  float lattice_velocity_max, lattice_stress_max;
  float pressure_min, pressure_max;
  float lattice_density_min, lattice_density_max;
  
  int pixels_x, pixels_y;
  
  vis_ctr_x     += steer_par[ 0 ];
  vis_ctr_y     += steer_par[ 1 ];
  vis_ctr_z     += steer_par[ 2 ];

  longitude      = steer_par[ 3 ];
  latitude       = steer_par[ 4 ];

  zoom           = steer_par[ 5 ];

  vis_brightness = steer_par[ 6 ];

  // The minimum value here is by default 0.0 all the time
  velocity_max   = steer_par[ 7 ];

  // The minimum value here is by default 0.0 all the time
  stress_max     = steer_par[ 8 ];

  pressure_min     = 80.F; // steer_par[ 9 ];
  pressure_max     = 110.F; // steer_par[ 10 ];

  vis_glyph_length = 1.F; // steer_par[ 11 ];

  pixels_x         = 512; // steer_par[ 12 ]; 
  pixels_y         = 512; // steer_par[ 13 ]; 
  
  vis_mouse_x      = (int)steer_par[ 14 ];
  vis_mouse_y      = (int)steer_par[ 15 ];

  lbm_terminate_simulation = (int)steer_par[ 16 ];

  *vis_perform_rendering   = (int)steer_par[ 17 ];

  visUpdateImageSize (pixels_x, pixels_y);
  
  visConvertThresholds (velocity_max, stress_max,
			pressure_min, pressure_max,
			&lattice_velocity_max, &lattice_stress_max,
			&lattice_density_min, &lattice_density_max,lbm);
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
  		 pixels_x, pixels_y,
  		 vis_ctr_x, vis_ctr_y, vis_ctr_z,
  		 5.F * vis->system_size,
  		 longitude, latitude,
  		 0.5F * (5.F * vis->system_size),
  		 zoom);
  
  vis_velocity_threshold_max_inv   = 1.F / lattice_velocity_max;
  vis_stress_threshold_max_inv     = 1.F / lattice_stress_max;
  vis_density_threshold_min        = 1.F / lattice_density_min;
  vis_density_threshold_minmax_inv = 1.F / (lattice_density_max - lattice_density_min);
}

