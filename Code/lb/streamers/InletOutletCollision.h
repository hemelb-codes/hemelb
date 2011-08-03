#ifndef HEMELB_LB_STREAMERS_INLETOUTLETCOLLISION_H
#define HEMELB_LB_STREAMERS_INLETOUTLETCOLLISION_H

#include "lb/streamers/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      class InletOutletCollision : public Collision
      {
        public:
          InletOutletCollision(distribn_t* iOutletDensityArray);

          ~InletOutletCollision();

          distribn_t getBoundaryDensityArray(const int index);

          virtual void AcceptCollisionVisitor(collisions::CollisionVisitor* v,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);

        private:
          distribn_t* mBoundaryDensityArray;
      };

    }
  }
}
#endif /* HEMELB_LB_STREAMERS_INLETOUTLETCOLLISION_H */
