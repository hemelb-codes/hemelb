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

          void PostStep(const bool iDoRayTracing,
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

          template<bool tDoRayTracing>
          void PostStepInternal(const double iOmega,
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

#endif /* HEMELB_LB_COLLISIONS_IMPLFINTERPOLATION_H */
