#ifndef HEMELB_LB_COLLISIONS_WALLCOLLISION_H
#define HEMELB_LB_COLLISIONS_WALLCOLLISION_H

#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class WallCollision : public Collision
      {
        protected:
          WallCollision();
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_WALLCOLLISION_H */
