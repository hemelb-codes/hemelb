#ifndef HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H

#include "lb/collisions/MidFluidCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplSimpleCollideAndStream : public MidFluidCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters &iLbmParams,
                            MinsAndMaxes* bMinimaAndMaxima,
                            LocalLatticeData &bLocalLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters &iLbmParams,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    LocalLatticeData &bLocalLatDat,
                                    hemelb::vis::Control *iControl);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H */
