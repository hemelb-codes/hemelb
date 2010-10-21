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
                                    const int iFirstIndex,
                                    const int iSiteCount,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    const Net* net,
                                    const double iStressType,
                                    const double iStressParam,
                                    hemelb::vis::Control *iControl);

          double* mBoundaryDensityArray;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H */
