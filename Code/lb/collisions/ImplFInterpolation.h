#ifndef HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H
#define HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplFInterpolation : public WallCollision
      {
        public:
          void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net);

          bool PostStep(double omega, int i, double *density, double *v_x, double *v_y,
            double *v_z, double f_neq[], Net* net);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H */
