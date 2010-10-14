#ifndef HEMELB_LB_COLLISIONS_COLLISION_H
#define HEMELB_LB_COLLISIONS_COLLISION_H

#include "net.h"
#include "D3Q15.h"

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
      };

    }
  }
}

#endif //HEMELB_LB_COLLISIONS_COLLISION_H
