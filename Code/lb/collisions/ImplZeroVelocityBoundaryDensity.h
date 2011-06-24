#ifndef HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/InletOutletWallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplZeroVelocityBoundaryDensity : public InletOutletWallCollision
      {
        public:
          ImplZeroVelocityBoundaryDensity(distribn_t* iBoundaryDensityArray);

          void DoCollisions(const bool iDoRayTracing,
                            const bool iDoEntropic,
                            const site_t iFirstIndex,
                            const site_t iSiteCount,
                            const LbmParameters* iLbmParams,
                            geometry::LatticeData* bLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const site_t iFirstIndex,
                                    const site_t iSiteCount,
                                    const LbmParameters* iLbmParams,
                                    geometry::LatticeData* bLatDat,
                                    hemelb::vis::Control *iControl);

          distribn_t* mBoundaryDensityArray;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H */
