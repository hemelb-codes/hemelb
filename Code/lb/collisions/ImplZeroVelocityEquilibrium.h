#ifndef HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplZeroVelocityEquilibrium : public WallCollision
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
                            const double iStressParam);

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
                                    const double iStressParam);
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H */
