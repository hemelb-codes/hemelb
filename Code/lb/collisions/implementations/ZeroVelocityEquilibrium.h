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

        template<bool tDoEntropic>
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void ZeroVelocityEquilibrium<tDoEntropic>::DoStreamAndCollide(WallCollision* mWallCollision,
                                                                      const site_t iFirstIndex,
                                                                      const site_t iSiteCount,
                                                                      const LbmParameters* iLbmParams,
                                                                      geometry::LatticeData* bLatDat,
                                                                      hemelb::vis::Control *iControl)
        {
          for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
          {
            distribn_t* lFOld = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
            distribn_t lFNeq[D3Q15::NUMVECTORS];
            distribn_t lDensity;

            lDensity = 0.0;

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lDensity += lFOld[ii];
            }

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lFNeq[ii] = lFOld[ii];
            }

            // Temporarily store FEq in lFNeq
            if (tDoEntropic)
              D3Q15::CalculateEntropicFeq(lDensity, 0.0, 0.0, 0.0, lFOld);
            else
              D3Q15::CalculateFeq(lDensity, 0.0, 0.0, 0.0, lFOld);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = lFOld[ii];
              lFNeq[ii] -= lFOld[ii];
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void ZeroVelocityEquilibrium<tDoEntropic>::DoPostStep(WallCollision* mWallCollision,
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
