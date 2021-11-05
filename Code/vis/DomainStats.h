// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_VIS_DOMAINSTATS_H
#define HEMELB_VIS_DOMAINSTATS_H

#include "constants.h"

namespace hemelb
{
  namespace vis
  {
    struct DomainStats
    {
      public:
        distribn_t velocity_threshold_max_inv;
        distribn_t stress_threshold_max_inv;
        distribn_t density_threshold_min;
        distribn_t density_threshold_minmax_inv;
        float physical_pressure_threshold_min;
        float physical_pressure_threshold_max;
        float physical_velocity_threshold_max;
        float physical_stress_threshold_max;
    };
  }
}

#endif /* HEMELB_VIS_DOMAINSTATS_H */

