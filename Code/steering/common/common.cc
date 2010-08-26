#include "constants.h"
#include "lb.h"

#include "vis/Control.h"
#include "vis/GlyphDrawer.h"

#include "steering/common/common.h"

namespace hemelb
{
  namespace steering
  {
    
    float steer_par[ STEERABLE_PARAMETERS + 1 ] = {
      0.0, 0.0, 0.0,    // scene center (dx,dy,dz)
      45.0, 45.0,          // longitude and latitude
      1.0, 0.03,        // zoom and brightness
      0.1, 0.1,         // velocity and stress ranges
      80.0, 120.0,      // Minimum pressure and maximum pressure for Colour mapping
      1.0,              // Glyph length
      512, 512,         // Rendered frame size, pixel x and pixel y
      -1.0, -1.0,         // x-y position of the mouse of the client
      0.0,              // signal useful to terminate the simulation
      0.0, 	         // Vis_mode 
      5.0,	         // vis_streaklines_per_pulsatile_period
      100.0,	         // vis_streakline_length
      0.0};             // doRendering

    void UpdateSteerableParameters(int *perform_rendering,
				   heme::vis::Control *visControl, LBM* lbm)
    {
      steer_par[ STEERABLE_PARAMETERS ] = *perform_rendering;

#ifndef NOMPI
      MPI_Bcast (steer_par, STEERABLE_PARAMETERS+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
      float longitude, latitude;
      float zoom;
      float lattice_velocity_max, lattice_stress_max;
      float lattice_density_min, lattice_density_max;
  
      int pixels_x, pixels_y;
  
      visControl->ctr_x     += steer_par[ 0 ];
      visControl->ctr_y     += steer_par[ 1 ];
      visControl->ctr_z     += steer_par[ 2 ];

      longitude      = steer_par[ 3 ];
      latitude       = steer_par[ 4 ];
  
      zoom           = steer_par[ 5 ];
  
      visControl->brightness = steer_par[ 6 ];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_velocity_threshold_max = steer_par[ 7 ];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_stress_threshold_max = steer_par[ 8 ];
  
      visControl->physical_pressure_threshold_min = steer_par[ 9 ];
      visControl->physical_pressure_threshold_max = steer_par[ 10 ];
  
      heme::vis::GlyphDrawer::glyph_length = steer_par[ 11 ];

      pixels_x         = steer_par[ 12 ]; 
      pixels_y         = steer_par[ 13 ]; 
  
      visControl->mouse_x      = int(steer_par[ 14 ]);
      visControl->mouse_y      = int(steer_par[ 15 ]);

      lbm_terminate_simulation = int(steer_par[ 16 ]);
  
      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      visControl->mode = int(steer_par[ 17 ]);
  
      visControl->streaklines_per_pulsatile_period = steer_par[ 18 ];
      visControl->streakline_length = steer_par[ 19 ];
  
      *perform_rendering = int(steer_par[ 20 ]);
  
      visControl->updateImageSize (pixels_x, pixels_y);
  
      lattice_density_min  = 
	lbm->lbmConvertPressureToLatticeUnits(visControl->physical_pressure_threshold_min) / Cs2;
      lattice_density_max  = lbm->lbmConvertPressureToLatticeUnits (visControl->physical_pressure_threshold_max) / Cs2;
      lattice_velocity_max = lbm->lbmConvertVelocityToLatticeUnits (visControl->physical_velocity_threshold_max);
      lattice_stress_max   = lbm->lbmConvertStressToLatticeUnits (visControl->physical_stress_threshold_max);  
  
      visControl->setProjection(pixels_x, pixels_y,
				longitude, latitude,
				zoom);
      
      visControl->density_threshold_min        = lattice_density_min;
      visControl->density_threshold_minmax_inv = 1.0F / (lattice_density_max - lattice_density_min);
      visControl->velocity_threshold_max_inv   = 1.0F / lattice_velocity_max;
      visControl->stress_threshold_max_inv     = 1.0F / lattice_stress_max;
    }

  }
}
