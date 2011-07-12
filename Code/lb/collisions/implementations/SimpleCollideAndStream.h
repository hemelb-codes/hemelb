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

        template<typename tCollisionOperator>
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

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void SimpleCollideAndStream<tCollisionOperator>::DoStreamAndCollide(MidFluidCollision* mMidFluidCollision,
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
            site_t siteIndex = iIndex - iFirstIndex;

            // Temporarily store f_eq in f_neq (rectified in next statement)
            tCollisionOperator::getSiteValues(lFOld, lDensity, lVx, lVy, lVz, lFNeq, siteIndex);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              // This also rectifies the lFNeq to actually store lFNeq and not lFEq
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                  += tCollisionOperator::getOperatorElement(lFOld[ii], lFNeq[ii], iLbmParams);
            }

            // lFOld is the post-collision, pre-streaming distribution
            tCollisionOperator::doPostCalculations(lFOld, bLatDat, siteIndex);

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

        template<typename tCollisionOperator>
        template<bool tDoRayTracing>
        void SimpleCollideAndStream<tCollisionOperator>::DoPostStep(MidFluidCollision* mMidFluidCollision,
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
