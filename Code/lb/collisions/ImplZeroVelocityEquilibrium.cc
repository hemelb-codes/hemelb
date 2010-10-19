#include "lb/collisions/ImplZeroVelocityEquilibrium.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      void ImplZeroVelocityEquilibrium::DoCollisions(const bool iDoRayTracing,
                                                     const double iOmega,
                                                     double iFOldAll[],
                                                     double iFNewAll[],
                                                     const int iFIdAll[],
                                                     const int iFirstIndex,
                                                     const int iSiteCount,
                                                     MinsAndMaxes* bMinimaAndMaxima,
                                                     const Net* net,
                                                     const double iStressType,
                                                     const double iStressParam)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                      iFirstIndex, iSiteCount,
                                      bMinimaAndMaxima, net, iStressType,
                                      iStressParam);
        }
        else
        {
          DoCollisionsInternal<false> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                       iFirstIndex, iSiteCount,
                                       bMinimaAndMaxima, net, iStressType,
                                       iStressParam);
        }
      }

      template<bool tDoRayTracing>
      void ImplZeroVelocityEquilibrium::DoCollisionsInternal(const double iOmega,
                                                             double iFOldAll[],
                                                             double iFNewAll[],
                                                             const int iFIdAll[],
                                                             const int iFirstIndex,
                                                             const int iSiteCount,
                                                             MinsAndMaxes* bMinimaAndMaxima,
                                                             const Net* net,
                                                             const double iStressType,
                                                             const double iStressParam)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = &iFOldAll[lIndex * D3Q15::NUMVECTORS];
          double lFNeq[15];
          double lDensity;

          lDensity = 0.0;

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lDensity += lFOld[ii];
          }

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq
          D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            iFNewAll[iFIdAll[lIndex * D3Q15::NUMVECTORS + ii]] = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing>(0.0, 0.0, 0.0, lIndex,
                                                       lFNeq, lDensity,
                                                       bMinimaAndMaxima, net,
                                                       iStressType,
                                                       iStressParam);
        }
      }
    }
  }
}
