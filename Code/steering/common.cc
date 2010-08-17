
#include "steering/steering.h"
#include "vis/rt.h"
#include "lb.h"

namespace steering {
  extern float steer_par[];
}

void UpdateSteerableParameters (int *perform_rendering,
				vis::Vis *vis, LBM* lbm)
{
  steering::steer_par[ STEERABLE_PARAMETERS ] = *perform_rendering;

#ifndef NOMPI
  MPI_Bcast (steering::steer_par, STEERABLE_PARAMETERS+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  float longitude, latitude;
  float zoom;
  float lattice_velocity_max, lattice_stress_max;
  float lattice_density_min, lattice_density_max;
  
  int pixels_x, pixels_y;
  
  vis::ctr_x     += steering::steer_par[ 0 ];
  vis::ctr_y     += steering::steer_par[ 1 ];
  vis::ctr_z     += steering::steer_par[ 2 ];

  longitude      = steering::steer_par[ 3 ];
  latitude       = steering::steer_par[ 4 ];
  
  zoom           = steering::steer_par[ 5 ];
  
  vis::brightness = steering::steer_par[ 6 ];

  // The minimum value here is by default 0.0 all the time
  vis::physical_velocity_threshold_max = steering::steer_par[ 7 ];

  // The minimum value here is by default 0.0 all the time
  vis::physical_stress_threshold_max = steering::steer_par[ 8 ];
  
  vis::physical_pressure_threshold_min = steering::steer_par[ 9 ];
  vis::physical_pressure_threshold_max = steering::steer_par[ 10 ];
  
  vis::GlyphDrawer::glyph_length = steering::steer_par[ 11 ];

  pixels_x         = steering::steer_par[ 12 ]; 
  pixels_y         = steering::steer_par[ 13 ]; 
  
  vis::mouse_x      = int(steering::steer_par[ 14 ]);
  vis::mouse_y      = int(steering::steer_par[ 15 ]);

  lbm_terminate_simulation = int(steering::steer_par[ 16 ]);
  
  // To swap between glyphs and streak line rendering...
  // 0 - Only display the isosurfaces (wall pressure and stress)
  // 1 - Isosurface and glyphs
  // 2 - Wall pattern streak lines
  vis::mode = int(steering::steer_par[ 17 ]);
  
  vis::streaklines_per_pulsatile_period = steering::steer_par[ 18 ];
  vis::streakline_length = steering::steer_par[ 19 ];
  
  *perform_rendering = int(steering::steer_par[ 20 ]);
  
  vis::visUpdateImageSize (pixels_x, pixels_y);
  
  lattice_density_min  = lbm->lbmConvertPressureToLatticeUnits (vis::physical_pressure_threshold_min) / Cs2;
  lattice_density_max  = lbm->lbmConvertPressureToLatticeUnits (vis::physical_pressure_threshold_max) / Cs2;
  lattice_velocity_max = lbm->lbmConvertVelocityToLatticeUnits (vis::physical_velocity_threshold_max);
  lattice_stress_max   = lbm->lbmConvertStressToLatticeUnits (vis::physical_stress_threshold_max);  
  
  vis::visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		      pixels_x, pixels_y,
		      vis::ctr_x, vis::ctr_y, vis::ctr_z,
		      5.F * vis->system_size,
		      longitude, latitude,
		      0.5F * (5.F * vis->system_size),
		      zoom);
  
  vis::density_threshold_min        = lattice_density_min;
  vis::density_threshold_minmax_inv = 1.0F / (lattice_density_max - lattice_density_min);
  vis::velocity_threshold_max_inv   = 1.0F / lattice_velocity_max;
  vis::stress_threshold_max_inv     = 1.0F / lattice_stress_max;
}

