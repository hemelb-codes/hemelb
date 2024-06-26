// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

#include <limits>
#include <cmath>
#include "units.h"

namespace hemelb
{

  constexpr unsigned int COLLISION_TYPES = 6;
  // TODO: when we hit C++20, use std::pi_v
  constexpr double PI = 3.14159265358979323846264338327950288;

  constexpr double mmHg_TO_PASCAL = 133.3223874;
  constexpr double DEFAULT_FLUID_DENSITY_Kg_per_m3 = 1000.0;
  constexpr double DEFAULT_FLUID_VISCOSITY_Pas = 0.004;

  /* This is the number of boundary types. It was 4, but the
   * "CHARACTERISTIC_BOUNDARY" type is never used and I don't know what it is
   * meant to be. It is also not used in the setup tool, so we will drop it,
   * setting BOUNDARIES to 3
   */
  constexpr sitedata_t BOUNDARIES = 3U;
  constexpr sitedata_t INLET_BOUNDARY = 0U;
  constexpr sitedata_t OUTLET_BOUNDARY = 1U;
  constexpr sitedata_t WALL_BOUNDARY = 2U;
  // const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

  constexpr unsigned int FLUID = 1U;
  constexpr unsigned int INLET = 2U;
  constexpr unsigned int OUTLET = 4U;
  constexpr unsigned int WALL = 8U;

  // square of the speed of sound
  constexpr double Cs2 = 1.0 / 3.0;
  // speed of the sound (sqrt is not constexpr)
  // TODO: when we hit C++20, use std::inv_sqrt3_v
  const double Cs = 1.0 / std::sqrt(3.0);

  // TODO almost certainly filth.
  constexpr distribn_t NO_VALUE = std::numeric_limits<distribn_t>::max();
  constexpr int SITE_OR_BLOCK_SOLID = std::numeric_limits<int>::min();
  constexpr int UNKNOWN_PROCESS = SITE_OR_BLOCK_SOLID + 1;

}

#endif //HEMELB_CONSTANTS_H
