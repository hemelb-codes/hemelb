#include "lb/collisions/ImplSimpleCollideAndStream.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void ImplSimpleCollideAndStream::DoCollisions(const bool iDoRayTracing,
                                                    const bool iDoEntropic,
                                                    const site_t iFirstIndex,
                                                    const site_t iSiteCount,
                                                    const LbmParameters *iLbmParams,
                                                    geometry::LatticeData *bLatDat,
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
      void ImplSimpleCollideAndStream::DoCollisionsInternal(const site_t iFirstIndex,
                                                            const site_t iSiteCount,
                                                            const LbmParameters *iLbmParams,
                                                            geometry::LatticeData *bLatDat,
                                                            hemelb::vis::Control *iControl)
      {
        for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
        {
          distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
          distribn_t lDensity, lVx, lVy, lVz;
          distribn_t lFNeq[D3Q15::NUMVECTORS];

          // Temporarily store f_eq in f_neq (rectified in next statement)
          D3Q15::CalculateDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz, lFNeq);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                += iLbmParams->Omega * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
          }

          UpdateMinsAndMaxes<tDoRayTracing> (lVx,
                                             lVy,
                                             lVz,
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

