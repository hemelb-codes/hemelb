#ifndef HEMELB_LB_COLLISIONS_COLLISION_H
#define HEMELB_LB_COLLISIONS_COLLISION_H

#include "net.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class Collision
      {
        public:
          virtual void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net) = 0;

        protected:
          void DensityAndVelocity(double f[], double *density, double *v_x, double *v_y,
            double *v_z);
          void CalculateFeq(double density, double v_x, double v_y, double v_z,
            double f_eq[]);
      };

    }
  }
}

#endif //HEMELB_LB_COLLISIONS_COLLISION_H
