#include "lb/collisions/ImplSimpleBounceBack.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleBounceBack::DoCollisions(const bool iDoRayTracing,
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
      void ImplSimpleBounceBack::DoCollisionsInternal(const double iOmega,
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
          double lFNeq[D3Q15::NUMVECTORS];
          double lVx, lVy, lVz, lDensity;

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          D3Q15::CalculateDensityVelocityFEq(lFNeq, lDensity, lVx, lVy, lVz,
                                             lFOld);

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            // The actual bounce-back lines, including streaming and collision. Basically swap the non-equilibrium components of f in each of the opposing pairs of directions.
            int lStreamTo = (net->HasBoundary(lIndex, ii))
              ? ( (lIndex * D3Q15::NUMVECTORS) + D3Q15::INVERSEDIRECTIONS[ii])
              : iFIdAll[lIndex * D3Q15::NUMVECTORS + ii];

            // Remember, oFNeq currently hold the equilibrium distribution. We
            // simultaneously use this and correct it, here.
            iFNewAll[lStreamTo] = lFOld[ii] += iOmega * (lFNeq[ii] -= lFOld[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, lIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima, net,
                                             iStressType, iStressParam);
        }
      }

    }
  }
}
