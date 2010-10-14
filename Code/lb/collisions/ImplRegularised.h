#ifndef HEMELB_LB_COLLISIONS_IMPLREGULARISED_H
#define HEMELB_LB_COLLISIONS_IMPLREGULARISED_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      // Implementation of BC3 from Chopard 2008
      class ImplRegularised : public WallCollision
      {
        public:
          void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net) = 0;
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLREGULARISED_H */
