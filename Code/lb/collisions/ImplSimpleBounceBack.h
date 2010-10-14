#ifndef HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H
#define HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplSimpleBounceBack : public WallCollision
      {
        public:
          void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net);
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H */
