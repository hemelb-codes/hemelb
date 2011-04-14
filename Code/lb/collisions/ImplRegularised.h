#ifndef HEMELB_LB_COLLISIONS_IMPLREGULARISED_H
#define HEMELB_LB_COLLISIONS_IMPLREGULARISED_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      // Implementation of BC3 from Chopard 2008
      class ImplRegularised : public WallCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters &iLbmParams,
                            geometry::LatticeData &bLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters &iLbmParams,
                                    geometry::LatticeData &bLatDat,
                                    hemelb::vis::Control *iControl);

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLREGULARISED_H */
