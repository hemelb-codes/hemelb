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
  float lattice_velocity_max, lattice_stress_max;
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
  vis_physical_velocity_threshold_max = steer_par[ 7 ];

  // The minimum value here is by default 0.0 all the time
  vis_physical_stress_threshold_max = steer_par[ 8 ];
  
  vis_physical_pressure_threshold_min = steer_par[ 9 ];
  vis_physical_pressure_threshold_max = steer_par[ 10 ];
  
  glyphDrawer::vis_glyph_length = steer_par[ 11 ];

  pixels_x         = steer_par[ 12 ]; 
  pixels_y         = steer_par[ 13 ]; 
  
  vis_mouse_x      = (int)steer_par[ 14 ];
  vis_mouse_y      = (int)steer_par[ 15 ];

  lbm_terminate_simulation = (int)steer_par[ 16 ];
  
  // To swap between glyphs and streak line rendering...
  // 0 - Only display the isosurfaces (wall pressure and stress)
  // 1 - Isosurface and glyphs
  // 2 - Wall pattern streak lines
  vis_mode = (int)steer_par[ 17 ];
  
  vis_streaklines_per_pulsatile_period = steer_par[ 18 ];
  vis_streakline_length = steer_par[ 19 ];
  
  *vis_perform_rendering = (int)steer_par[ 20 ];
  
  visUpdateImageSize (pixels_x, pixels_y);
  
  lattice_density_min  = lbm->lbmConvertPressureToLatticeUnits (vis_physical_pressure_threshold_min) / Cs2;
  lattice_density_max  = lbm->lbmConvertPressureToLatticeUnits (vis_physical_pressure_threshold_max) / Cs2;
  lattice_velocity_max = lbm->lbmConvertVelocityToLatticeUnits (vis_physical_velocity_threshold_max);
  lattice_stress_max   = lbm->lbmConvertStressToLatticeUnits (vis_physical_stress_threshold_max);  
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
  		 pixels_x, pixels_y,
  		 vis_ctr_x, vis_ctr_y, vis_ctr_z,
  		 5.F * vis->system_size,
  		 longitude, latitude,
  		 0.5F * (5.F * vis->system_size),
  		 zoom);
  
  vis_density_threshold_min        = lattice_density_min;
  vis_density_threshold_minmax_inv = 1.0F / (lattice_density_max - lattice_density_min);
  vis_velocity_threshold_max_inv   = 1.0F / lattice_velocity_max;
  vis_stress_threshold_max_inv     = 1.0F / lattice_stress_max;
}

