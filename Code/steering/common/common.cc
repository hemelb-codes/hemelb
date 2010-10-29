#include "constants.h"

#include "steering/common/common.h"

namespace hemelb
{
  namespace steering
  {

    float steer_par[STEERABLE_PARAMETERS + 1] = { 0.0, 0.0, 0.0, // scene center (dx,dy,dz)
                                                  45.0, 45.0, // longitude and latitude
                                                  1.0, 0.03, // zoom and brightness
                                                  0.1, 0.1, // velocity and stress ranges
                                                  80.0, 120.0, // Minimum pressure and maximum pressure for Colour mapping
                                                  1.0, // Glyph length
                                                  512, 512, // Rendered frame size, pixel x and pixel y
                                                  -1.0, -1.0, // x-y position of the mouse of the client
                                                  0.0, // signal useful to terminate the simulation
                                                  0.0, // Vis_mode
                                                  5.0, // vis_streaklines_per_pulsatile_period
                                                  100.0, // vis_streakline_length
                                                  0.0 }; // doRendering


  }
}
