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
  
  
  vis_ctr_x     += steer_par[ 0 ];
  vis_ctr_y     += steer_par[ 1 ];
  vis_ctr_z     += steer_par[ 2 ];
  longitude      = steer_par[ 3 ];
  latitude       = steer_par[ 4 ];
  zoom           = steer_par[ 5 ];
  vis_brightness = steer_par[ 6 ];
  velocity_max   = steer_par[ 7 ];
  stress_max     = steer_par[ 8 ];
  
  vis_mouse_x              = (int)steer_par[  9 ];
  vis_mouse_y              = (int)steer_par[ 10 ];
  lbm_terminate_simulation = (int)steer_par[ 11 ];
  *vis_perform_rendering   = (int)steer_par[ 12 ];
  
  visConvertThresholds (velocity_max, stress_max,
			&lattice_velocity_max, &lattice_stress_max, lbm);
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
  		 PIXELS_X, PIXELS_Y,
  		 vis_ctr_x, vis_ctr_y, vis_ctr_z,
  		 5.F * vis->system_size,
  		 longitude, latitude,
  		 0.5F * (5.F * vis->system_size),
  		 zoom);
  
  vis_velocity_threshold_max_inv = 1.0/lattice_velocity_max;
  vis_stress_threshold_max_inv   = 1.0/lattice_stress_max;
}
