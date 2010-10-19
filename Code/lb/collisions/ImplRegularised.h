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
                            const double iOmega,
                            double iFOldAll[],
                            double iFNewAll[],
                            const int iFIdAll[],
                            const int iFirstIndex,
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
                                    const int iFirstIndex,
                                    const int iSiteCount,
                                    MinsAndMaxes* bMinimaAndMaxima,
                                    const Net* net,
                                    const double iStressType,
                                    const double iStressParam);

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLREGULARISED_H */
