#ifndef HEMELB_LB_COLLISIONS_COLLISION_H
#define HEMELB_LB_COLLISIONS_COLLISION_H

#include "vis/Control.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "lb/BoundaryComms.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class CollisionVisitor; // To deal with circular dependency of collision and visitor

      class Collision
      {
        public:
          virtual ~Collision();

          virtual void AcceptCollisionVisitor(CollisionVisitor* v,
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

#endif /* HEMELB_LB_COLLISIONS_COLLISION_H */
