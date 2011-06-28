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
        public:
          virtual void Accept(Visitor* v,
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
#endif /* HEMELB_LB_COLLISIONS_WALLCOLLISION_H */
