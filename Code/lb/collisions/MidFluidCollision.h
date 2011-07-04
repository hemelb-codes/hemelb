#ifndef HEMELB_LB_COLLISIONS_MIDFLUIDCOLLISION_H
#define HEMELB_LB_COLLISIONS_MIDFLUIDCOLLISION_H

#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class MidFluidCollision : public Collision
      {
        public:
          virtual ~MidFluidCollision();

          virtual void AcceptCollisionVisitor(CollisionVisitor* v,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_MIDFLUIDCOLLISION_H */
