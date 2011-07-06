#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYEQUILIBRIUM_H

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
        class ZeroVelocityEquilibrium : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(WallCollision* mWallCollision,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

            template<bool tDoRayTracing>
            static void DoPostStep(WallCollision* mWallCollision,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl);

        };

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void ZeroVelocityEquilibrium<tCollisionOperator>::DoStreamAndCollide(WallCollision* mWallCollision,
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
            distribn_t lDensity;

            lDensity = 0.0;

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lDensity += lFOld[ii];
            }

            tCollisionOperator::getBoundarySiteValues(lFOld, lDensity, 0.0, 0.0, 0.0, lFEq, lIndex - iFirstIndex);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFEq[ii];
              lFNeq[ii] = lFOld[ii] - lFEq[ii];
              lFOld[ii] = lFEq[ii];
            }

            UpdateMinsAndMaxes<tDoRayTracing> (0.0,
                                               0.0,
                                               0.0,
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
        void ZeroVelocityEquilibrium<tCollisionOperator>::DoPostStep(WallCollision* mWallCollision,
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

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_ZEROVELOCITYEQUILIBRIUM */
