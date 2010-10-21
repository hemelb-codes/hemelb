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
                            const double iOmega,
                            double iFOldAll[],
                            double iFNewAll[],
                            const int iFIdAll[],
                            int iFirstIndex,
                            const int iSiteCount,
                            MinsAndMaxes* bMinimaAndMaxima,
                            const Net* net,
                            const double iStressType,
                            const double iStressParam,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const double iOmega,
                                    double iFOldAll[],
                                    double iFNewAll[],
                                    const int iFIdAll[],
                                    int iFirstIndex,
                                    const int iSiteCount,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    const Net* net,
                                    const double iStressType,
                                    const double iStressParam,
                                    hemelb::vis::Control *iControl);
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLECOLLIDEANDSTREAM_H */
