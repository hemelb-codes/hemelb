#ifndef HEMELB_LB_STREAMERS_COLLISION_H
#define HEMELB_LB_STREAMERS_COLLISION_H

#include "vis/Control.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "lb/boundaries/BoundaryValues.h"


namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      class CollisionVisitor;
    }

    namespace streamers
    {
      class Collision
      {
        public:
          virtual ~Collision();

          virtual void AcceptCollisionVisitor(hemelb::lb::collisions::CollisionVisitor* v,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl) = 0;

      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_COLLISION_H */
