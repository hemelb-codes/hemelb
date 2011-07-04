#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_SIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_SIMPLECOLLIDEANDSTREAM_H

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
        class SimpleCollideAndStream : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(MidFluidCollision* mMidFluidCollision,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

            template<bool tDoRayTracing>
            static void DoPostStep(MidFluidCollision* mMidFluidCollision,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl);

        };

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void SimpleCollideAndStream<tDoEntropic>::DoStreamAndCollide(MidFluidCollision* mMidFluidCollision,
                                                                     const site_t iFirstIndex,
                                                                     const site_t iSiteCount,
                                                                     const LbmParameters* iLbmParams,
                                                                     geometry::LatticeData* bLatDat,
                                                                     hemelb::vis::Control *iControl)
        {
          for (site_t iIndex = iFirstIndex; iIndex < (iFirstIndex + iSiteCount); iIndex++)
          {
            distribn_t* lFOld = bLatDat->GetFOld(iIndex * D3Q15::NUMVECTORS);
            distribn_t lDensity, lVx, lVy, lVz;
            distribn_t lFNeq[D3Q15::NUMVECTORS];
            double alpha;

            // Temporarily store f_eq in f_neq (rectified in next statement)
            if (tDoEntropic)
            {
              D3Q15::CalculateEntropicDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz, lFNeq);
              alpha = getAlpha(lFOld, lFNeq); // lFNeq is actually lFeq at the moment
            }
            else
            {
              D3Q15::CalculateDensityVelocityFEq(lFOld, lDensity, lVx, lVy, lVz, lFNeq);
            }

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              if (tDoEntropic)
              {
                * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii] += alpha
                    * iLbmParams->Beta * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
              }
              else
              {
                * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                    += iLbmParams->Omega * (lFNeq[ii] = lFOld[ii] - lFNeq[ii]);
              }
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

        template<bool tDoEntropic>
        template<bool tDoRayTracing>
        void SimpleCollideAndStream<tDoEntropic>::DoPostStep(MidFluidCollision* mMidFluidCollision,
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

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_SIMPLECOLLIDEANDSTREAM */
