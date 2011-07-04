#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYBOUNDARYDENSITY_H

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
        class ZeroVelocityBoundaryDensity : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(InletOutletWallCollision* mInletOutletWallCollision,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

            template<bool tDoRayTracing>
            static void DoPostStep(InletOutletWallCollision* mInletOutletWallCollision,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl);

        };

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void ZeroVelocityBoundaryDensity<tDoEntropic>::DoStreamAndCollide(InletOutletWallCollision* mInletOutletWallCollision,
                                                                          const site_t iFirstIndex,
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

            lDensity
                = (*mInletOutletWallCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(iIndex));

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lFNeq[ii] = lFOld[ii];
            }

            // Temporarily store FEq in FNeq
            if (tDoEntropic)
              D3Q15::CalculateEntropicFeq(lDensity, 0.0, 0.0, 0.0, lFOld);
            else
              D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii];
              lFNeq[ii] -= lFOld[ii];
            }

            UpdateMinsAndMaxes<tDoRayTracing> (0.0,
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void ZeroVelocityBoundaryDensity<tDoEntropic>::DoPostStep(InletOutletWallCollision* mInletOutletWallCollision,
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

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYBOUNDARYDENSITY */
