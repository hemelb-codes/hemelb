#include "lb/collisions/ImplSimpleCollideAndStream.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleCollideAndStream::DoCollisions(const bool iDoRayTracing,
                                                    const double iOmega,
                                                    double iFOldAll[],
                                                    double iFNewAll[],
                                                    const int iFIdAll[],
                                                    const int iFirstIndex,
                                                    const int iSiteCount,
                                                    MinsAndMaxes* bMinimaAndMaxima,
                                                    const Net* net,
                                                    const double iStressType,
                                                    const double iStressParam,
                                                    hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                      iFirstIndex, iSiteCount,
                                      bMinimaAndMaxima, net, iStressType,
                                      iStressParam, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iOmega, iFOldAll, iFNewAll, iFIdAll,
                                       iFirstIndex, iSiteCount,
                                       bMinimaAndMaxima, net, iStressType,
                                       iStressParam, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplSimpleCollideAndStream::DoCollisionsInternal(const double iOmega,
                                                            double iFOldAll[],
                                                            double iFNewAll[],
                                                            const int iFIdAll[],
                                                            const int iFirstIndex,
                                                            const int iSiteCount,
                                                            MinsAndMaxes* bMinimaAndMaxima,
                                                            const Net* net,
                                                            const double iStressType,
                                                            const double iStressParam,
                                                            hemelb::vis::Control *iControl)
      {
        for (int iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          double *lFOld = &iFOldAll[iIndex * D3Q15::NUMVECTORS];
          double lDensity, lVx, lVy, lVz;
          double lFNeq[15];

          // Temporarily store f_eq in f_neq (rectified in next statement)
          D3Q15::CalculateDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz,
                                             lFNeq);

          for (int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            iFNewAll[iFIdAll[iIndex * D3Q15::NUMVECTORS + ii]] = lFOld[ii]
                += iOmega * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, iIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima, net,
                                             iStressType, iStressParam, iControl);
        }
      }

    }
  }
}

