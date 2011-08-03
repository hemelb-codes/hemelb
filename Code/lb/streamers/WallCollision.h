#ifndef HEMELB_LB_STREAMERS_WALLCOLLISION_H
#define HEMELB_LB_STREAMERS_WALLCOLLISION_H

#include "lb/streamers/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      class WallCollision : public Collision
      {
        public:
          virtual void AcceptCollisionVisitor(collisions::CollisionVisitor* v,
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
#endif /* HEMELB_LB_STREAMERS_WALLCOLLISION_H */
