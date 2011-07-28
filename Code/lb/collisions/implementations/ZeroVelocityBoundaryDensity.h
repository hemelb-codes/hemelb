#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/implementations/Implementation.h"
#include "lb/collisions/implementations/CollisionOperator.h"

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
                                           tCollisionOperator* iCollisionOperator,
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
                                                                                 tCollisionOperator* iCollisionOperator,
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
            site_t siteIndex = iIndex - iFirstIndex;

            lDensity
                = (*mInletOutletWallCollision).getBoundaryDensityArray(bLatDat->GetBoundaryId(iIndex));

            iCollisionOperator->getBoundarySiteValues(lFOld,
                                                      lDensity,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      lFEq,
                                                      siteIndex);

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
