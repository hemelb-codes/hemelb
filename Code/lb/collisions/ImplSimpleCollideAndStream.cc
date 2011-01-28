#include "lb/collisions/ImplSimpleCollideAndStream.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleCollideAndStream::DoCollisions(const bool iDoRayTracing,
                                                    const int iFirstIndex,
                                                    const int iSiteCount,
                                                    const LbmParameters &iLbmParams,
                                                    MinsAndMaxes &bMinimaAndMaxima,
                                                    LocalLatticeData &bLocalLatDat,
                                                    hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams,
                                      bMinimaAndMaxima, bLocalLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams,
                                       bMinimaAndMaxima, bLocalLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplSimpleCollideAndStream::DoCollisionsInternal(const int iFirstIndex,
                                                            const int iSiteCount,
                                                            const LbmParameters &iLbmParams,
                                                            MinsAndMaxes &bMinimaAndMaxima,
                                                            LocalLatticeData &bLocalLatDat,
                                                            hemelb::vis::Control *iControl)
      {
        for (int iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          double *lFOld = &bLocalLatDat .FOld[iIndex * D3Q15::NUMVECTORS];
          double lDensity, lVx, lVy, lVz;
          double lFNeq[D3Q15::NUMVECTORS];

          // Temporarily store f_eq in f_neq (rectified in next statement)
          D3Q15::CalculateDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz,
                                             lFNeq);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            bLocalLatDat .FNew[bLocalLatDat.GetStreamedIndex(iIndex, ii)]
                = lFOld[ii] += iLbmParams.Omega * (lFNeq[ii] = lFOld[ii]
                    - lFNeq[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, iIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima,
                                             bLocalLatDat, iLbmParams, iControl);
        }
      }

    }
  }
}

