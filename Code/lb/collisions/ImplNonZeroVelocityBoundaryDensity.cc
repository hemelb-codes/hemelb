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
      void ImplNonZeroVelocityBoundaryDensity::DoCollisionsInternal(const int iFirstIndex,
                                                                    const int iSiteCount,
                                                                    const LbmParameters &iLbmParams,
                                                                    MinsAndMaxes &bMinimaAndMaxima,
                                                                    LocalLatticeData &bLocalLatDat,
                                                                    hemelb::vis::Control *iControl)
      {
        for (int lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          double *lFOld = &bLocalLatDat.FOld[lIndex * D3Q15::NUMVECTORS];
          double lFNeq[15];
          double lVx, lVy, lVz, lDummyDensity, lDensity;

          lDensity = mBoundaryDensityArray[bLocalLatDat.GetBoundaryId(lIndex)];

          D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy,
                                             lVz);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq (rectified later).
          D3Q15::CalculateFeq(lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            bLocalLatDat.FNew[bLocalLatDat.GetStreamedIndex(lIndex, ii)]
                = lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx, lVy, lVz, lIndex, lFNeq,
                                             lDensity, bMinimaAndMaxima,
                                             bLocalLatDat, iLbmParams, iControl);
        }
      }

    }
  }
}
