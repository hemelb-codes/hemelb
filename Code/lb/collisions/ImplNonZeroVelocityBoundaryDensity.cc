#include "lb/collisions/ImplNonZeroVelocityBoundaryDensity.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      ImplNonZeroVelocityBoundaryDensity::ImplNonZeroVelocityBoundaryDensity(double* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void ImplNonZeroVelocityBoundaryDensity::DoCollisions(const bool iDoRayTracing,
                                                            const int iFirstIndex,
                                                            const int iSiteCount,
                                                            const LbmParameters &iLbmParams,
                                                            MinsAndMaxes &bMinimaAndMaxima,
                                                            geometry::LatticeData &bLatDat,
                                                            hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                      bLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bMinimaAndMaxima,
                                       bLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplNonZeroVelocityBoundaryDensity::DoCollisionsInternal(const int iFirstIndex,
                                                                    const int iSiteCount,
                                                                    const LbmParameters &iLbmParams,
                                                                    MinsAndMaxes &bMinimaAndMaxima,
                                                                    geometry::LatticeData &bLatDat,
                                                                    hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = bLatDat.GetFOld(lIndex * D3Q15::NUMVECTORS);
          double lFNeq[15];
          double lVx, lVy, lVz, lDummyDensity, lDensity;

          lDensity = mBoundaryDensityArray[bLatDat.GetBoundaryId(lIndex)];

          D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy, lVz);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq (rectified later).
          D3Q15::CalculateFeq(lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat.GetFNew(bLatDat.GetStreamedIndex(lIndex, ii))) = lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, lIndex, lFNeq, lDensity,
                                             bMinimaAndMaxima, bLatDat, iLbmParams, iControl);
        }
      }

    }
  }
}
