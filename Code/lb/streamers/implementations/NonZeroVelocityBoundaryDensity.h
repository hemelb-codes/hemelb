#ifndef HEMELB_LB_STREAMERS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_STREAMERS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/streamers/Implementation.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      namespace implementations
      {

        template<typename tCollisionOperator>
        class NonZeroVelocityBoundaryDensity : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(InletOutletCollision* mInletOutletCollision,
                                           tCollisionOperator* iCollisionOperator,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

            template<bool tDoRayTracing>
            static void DoPostStep(InletOutletCollision* mInletOutletCollision,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl);

        };

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void NonZeroVelocityBoundaryDensity<tCollisionOperator>::DoStreamAndCollide(InletOutletCollision* mInletOutletCollision,
                                                                                    tCollisionOperator* iCollisionOperator,
                                                                                    const site_t iFirstIndex,
                                                                                    const site_t iSiteCount,
                                                                                    const LbmParameters* iLbmParams,
                                                                                    geometry::LatticeData* bLatDat,
                                                                                    hemelb::vis::Control *iControl)
        {
          for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
          {
            distribn_t* lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
            distribn_t lFNeq[D3Q15::NUMVECTORS], lFEq[D3Q15::NUMVECTORS];
            distribn_t lVx, lVy, lVz, lDummyDensity, lDensity;
            site_t siteIndex = lIndex - iFirstIndex;

            lDensity
                = mInletOutletCollision->getBoundaryDensityArray(bLatDat->GetBoundaryId(lIndex));

            D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy, lVz);

            iCollisionOperator->getBoundarySiteValues(lFOld,
                                                      lDensity,
                                                      lVx,
                                                      lVy,
                                                      lVz,
                                                      lFEq,
                                                      siteIndex);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFEq[ii];
              lFNeq[ii] = lFOld[ii] - lFEq[ii];
              lFOld[ii] = lFEq[ii];
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

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void NonZeroVelocityBoundaryDensity<tCollisionOperator>::DoPostStep(InletOutletCollision* mInletOutletCollision,
                                                                            const site_t iFirstIndex,
                                                                            const site_t iSiteCount,
                                                                            const LbmParameters* iLbmParams,
                                                                            geometry::LatticeData* bLatDat,
                                                                            hemelb::vis::Control *iControl)
        {

        }

      }
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H */
