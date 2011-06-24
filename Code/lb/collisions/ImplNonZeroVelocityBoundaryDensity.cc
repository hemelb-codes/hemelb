#include "lb/collisions/ImplNonZeroVelocityBoundaryDensity.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      ImplNonZeroVelocityBoundaryDensity::ImplNonZeroVelocityBoundaryDensity(distribn_t* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void ImplNonZeroVelocityBoundaryDensity::DoCollisions(const bool iDoRayTracing,
                                                            const bool iDoEntropic,
                                                            const site_t iFirstIndex,
                                                            const site_t iSiteCount,
                                                            const LbmParameters* iLbmParams,
                                                            geometry::LatticeData* bLatDat,
                                                            hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          DoCollisionsInternal<true> (iFirstIndex, iSiteCount, iLbmParams, bLatDat, iControl);
        }
        else
        {
          DoCollisionsInternal<false> (iFirstIndex, iSiteCount, iLbmParams, bLatDat, iControl);
        }
      }

      template<bool tDoRayTracing>
      void ImplNonZeroVelocityBoundaryDensity::DoCollisionsInternal(const site_t iFirstIndex,
                                                                    const site_t iSiteCount,
                                                                    const LbmParameters* iLbmParams,
                                                                    geometry::LatticeData* bLatDat,
                                                                    hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[15];
          distribn_t lVx, lVy, lVz, lDummyDensity, lDensity;

          lDensity = mBoundaryDensityArray[bLatDat->GetBoundaryId(lIndex)];

          D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy, lVz);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in lFNeq (rectified later).
          D3Q15::CalculateFeq(lDensity, lVx, lVy, lVz, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFOld[ii];
          }

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] -= lFOld[ii];
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx,
                                             lVy,
                                             lVz,
                                             lIndex,
                                             lFNeq,
                                             lDensity,
                                             bLatDat,
                                             iLbmParams,
                                             iControl);
        }
      }

    }
  }
}
