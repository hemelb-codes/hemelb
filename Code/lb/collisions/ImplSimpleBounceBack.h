#ifndef HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H
#define HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplSimpleBounceBack : public WallCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const bool iDoEntropic,
                            const site_t iFirstIndex,
                            const site_t iSiteCount,
                            const LbmParameters *iLbmParams,
                            geometry::LatticeData *bLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const site_t iFirstIndex,
                                    const site_t iSiteCount,
                                    const LbmParameters *iLbmParams,
                                    geometry::LatticeData *bLatDat,
                                    hemelb::vis::Control *iControl);

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H */
