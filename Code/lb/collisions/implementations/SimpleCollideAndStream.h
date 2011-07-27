#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_SIMPLECOLLIDEANDSTREAM_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_SIMPLECOLLIDEANDSTREAM_H

#include "lb/collisions/implementations/Implementation.h"
#include "lb/collisions/implementations/HFunction.h"

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
            distribn_t lFNeq[D3Q15::NUMVECTORS], lFEq[D3Q15::NUMVECTORS];
            double alpha;

            tCollisionOperator::getSiteValues(lFOld, lDensity, lVx, lVy, lVz, lFEq, iIndex
                - iFirstIndex);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              lFNeq[ii] = lFOld[ii] - lFEq[ii];
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(iIndex, ii))) = lFOld[ii]
                  += tCollisionOperator::getOperatorElement(lFOld[ii], lFNeq[ii], iLbmParams);
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
