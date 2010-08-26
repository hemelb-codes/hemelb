#ifndef HEMELB_STEERING_COMMON_H
#define HEMELB_STEERING_COMMON_H

#include "lb.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {
    void UpdateSteerableParameters(int *perform_rendering,
				   heme::vis::Control *visControl, LBM* lbm);

    extern float steer_par[ STEERABLE_PARAMETERS + 1 ];
  }
}
#endif // HEMELB_STEERING_COMMON_H
