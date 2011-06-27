#include "lb/collisions/ImplZeroVelocityBoundaryDensity.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      ImplZeroVelocityBoundaryDensity::ImplZeroVelocityBoundaryDensity(distribn_t* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void ImplZeroVelocityBoundaryDensity::DoCollisions(const bool iDoRayTracing,                                                         const site_t iFirstIndex,
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

      // Collision + streaming for fluid lattice sites and adjacent to the outlet and the wall.
      template<bool tDoRayTracing>
      void ImplZeroVelocityBoundaryDensity::DoCollisionsInternal(const site_t iFirstIndex,
                                                                 const site_t iSiteCount,
                                                                 const LbmParameters* iLbmParams,
                                                                 geometry::LatticeData* bLatDat,
                                                                 hemelb::vis::Control *iControl)
      {
        for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
          distribn_t lFNeq[D3Q15::NUMVECTORS];
          distribn_t lDensity;

          lDensity = mBoundaryDensityArray[bLatDat->GetBoundaryId(iIndex)];

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            lFNeq[ii] = lFOld[ii];
          }

          // Temporarily store FEq in FNeq
          D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii];
            lFNeq[ii] -= lFOld[ii];
          }

          Collision::UpdateMinsAndMaxes<tDoRayTracing>(0.0,
                                                       0.0,
                                                       0.0,
                                                       iIndex,
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
