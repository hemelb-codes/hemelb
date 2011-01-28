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
          ImplZeroVelocityBoundaryDensity(double* iBoundaryDensityArray);

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

#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYBOUNDARYDENSITY_H */
