#ifndef HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H
#define HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplFInterpolation : public WallCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters* iLbmParams,
                            geometry::LatticeData* bLatDat,
                            hemelb::vis::Control *iControl);

          void PostStep(const bool iDoRayTracing,
                        const int iFirstIndex,
                        const int iSiteCount,
                        const LbmParameters* iLbmParams,
                        geometry::LatticeData* bLatDat,
                        hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters* iLbmParams,
                                    geometry::LatticeData* bLatDat,
                                    hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void PostStepInternal(const int iFirstIndex,
                                const int iSiteCount,
                                const LbmParameters* iLbmParams,
                                geometry::LatticeData* bLatDat,
                                hemelb::vis::Control *iControl);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H */
