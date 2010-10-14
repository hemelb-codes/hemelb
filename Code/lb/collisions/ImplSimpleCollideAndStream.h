#ifndef HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H

#include "lb/collisions/MidFluidCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplSimpleCollideAndStream : public MidFluidCollision
      {
          void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H */
