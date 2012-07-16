// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_CONSTANTS_H
#define HEMELB_CONSTANTS_H

#include <limits>
#include "units.h"

namespace hemelb
{

  const unsigned int COLLISION_TYPES = 6;
  const double PI = 3.14159265358979323846264338327950288;
  const double DEG_TO_RAD = (PI / 180.0);
  // TODO this was used for a convergence test - we could reinstate that at some point.
  const double EPSILON = 1.0e-30;

  const double REFERENCE_PRESSURE_mmHg = 80.0;
  const double mmHg_TO_PASCAL = 133.3223874;
  const double BLOOD_DENSITY_Kg_per_m3 = 1000.0;
  const double BLOOD_VISCOSITY_Pa_s = 0.004;

  /* This is the number of boundary types. It was 4, but the
   * "CHARACTERISTIC_BOUNDARY" type is never used and I don't know what it is
   * meant to be. It is also not used in the setup tool, so we will drop it,
   * setting BOUNDARIES to 3
   */
  const sitedata_t BOUNDARIES = 3U;
  const sitedata_t INLET_BOUNDARY = 0U;
  const sitedata_t OUTLET_BOUNDARY = 1U;
  const sitedata_t WALL_BOUNDARY = 2U;
  // const unsigned int CHARACTERISTIC_BOUNDARY = 3U;

  const unsigned int FLUID = 1U;
  const unsigned int INLET = 2U;
  const unsigned int OUTLET = 4U;
  const unsigned int EDGE = 8U;

  // square of the speed of sound
  const double Cs2 = 1.0 / 3.0;

  // TODO almost certainly filth.
  const distribn_t NO_VALUE = std::numeric_limits<distribn_t>::max();
  const int BIG_NUMBER2 = 1 << 30;
}

#endif //HEMELB_CONSTANTS_H
