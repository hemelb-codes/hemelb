#ifndef HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/InletOutletCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplNonZeroVelocityBoundaryDensity : public InletOutletCollision
      {
        public:
          ImplNonZeroVelocityBoundaryDensity(double* iBounaryDensityArray);

          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters &iLbmParams,
                            MinsAndMaxes &bMinimaAndMaxima,
                            LocalLatticeData &bLocalLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters &iLbmParams,
                                    MinsAndMaxes &bMinimaAndMaxima,
                                    LocalLatticeData &bLocalLatDat,
                                    hemelb::vis::Control *iControl);

          double* mBoundaryDensityArray;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H */
