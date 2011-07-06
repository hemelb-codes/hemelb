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

        template<typename tCollisionOperator>
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

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void ZeroVelocityBoundaryDensity<tCollisionOperator>::DoStreamAndCollide(InletOutletWallCollision* mInletOutletWallCollision,
                                                                          const site_t iFirstIndex,
                                                                          const site_t iSiteCount,
                                                                          const LbmParameters* iLbmParams,
                                                                          geometry::LatticeData* bLatDat,
                                                                          hemelb::vis::Control *iControl)
        {
          for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
          {
            distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
            distribn_t lFNeq[D3Q15::NUMVECTORS], lFEq[D3Q15::NUMVECTORS];
            distribn_t lDensity;

            lDensity
                = (*mInletOutletWallCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(iIndex));

            tCollisionOperator::getBoundarySiteValues(lFOld, lDensity, 0.0, 0.0, 0.0, lFEq, iIndex - iFirstIndex);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFEq[ii];
              lFNeq[ii] = lFOld[ii] - lFEq[ii];
              lFOld[ii] = lFEq[ii];
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

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void ZeroVelocityBoundaryDensity<tCollisionOperator>::DoPostStep(InletOutletWallCollision* mInletOutletWallCollision,
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
