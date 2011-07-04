#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/implementations/Implementation.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        template<bool tDoEntropic>
        class NonZeroVelocityBoundaryDensity : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(InletOutletCollision* mInletOutletCollision,
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void NonZeroVelocityBoundaryDensity<tDoEntropic>::DoStreamAndCollide(InletOutletCollision* mInletOutletCollision,
                                                                             const site_t iFirstIndex,
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

            lDensity
                = (*mInletOutletCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(lIndex));

            D3Q15::CalculateDensityAndVelocity(lFOld, lDummyDensity, lVx, lVy, lVz);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lFNeq[ii] = lFOld[ii];
            }

            // Temporarily store FEq in lFNeq (rectified later).
            if (tDoEntropic)
              D3Q15::CalculateEntropicFeq(lDensity, lVx, lVy, lVz, lFOld);
            else
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void NonZeroVelocityBoundaryDensity<tDoEntropic>::DoPostStep(InletOutletCollision* mInletOutletCollision,
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

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_NONZEROVELOCITYBOUNDARYDENSITY_H */
