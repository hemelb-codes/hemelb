#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_FINTERPOLATION_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_FINTERPOLATION_H

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

        class FInterpolation : public Implementation
        {

          public:
            template<bool tDoRayTracing>
            static void DoStreamAndCollide(WallCollision* mWallCollision,
                                           CollisionOperator* iCollisionOperator,
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

        template<bool tDoRayTracing>
        void FInterpolation::DoStreamAndCollide(WallCollision* mWallCollision,
                                                CollisionOperator* iCollisionOperator,
                                                const site_t iFirstIndex,
                                                const site_t iSiteCount,
                                                const LbmParameters* iLbmParams,
                                                geometry::LatticeData* bLatDat,
                                                hemelb::vis::Control *iControl)
        {
          for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
          {
            distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
            distribn_t density, v_x, v_y, v_z, f_neq[15];
            // Temporarily store f_eq in f_neq. Rectified later.
            D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

            for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
            {
              * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = f[ii]
                  += iLbmParams->Omega * (f_neq[ii] = f[ii] - f_neq[ii]);
            }

            UpdateMinsAndMaxes<tDoRayTracing> (v_x,
                                               v_y,
                                               v_z,
                                               lIndex,
                                               f_neq,
                                               density,
                                               bLatDat,
                                               iLbmParams,
                                               iControl);
          }
        }

        template<bool tDoRayTracing>
        void FInterpolation::DoPostStep(WallCollision* mWallCollision,
                                        const site_t iFirstIndex,
                                        const site_t iSiteCount,
                                        const LbmParameters* iLbmParams,
                                        geometry::LatticeData* bLatDat,
                                        hemelb::vis::Control *iControl)
        {
          for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
          {
            // Handle odd indices, then evens - it's slightly easier to take the odd
            // and even cases separately.
            for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
            {
              if (bLatDat->HasBoundary(lIndex, l))
              {
                double twoQ = 2.0 * bLatDat->GetCutDistance(lIndex, l);

                * (bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l]))
                    = (twoQ < 1.0)
                      ? (*bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + l) + twoQ
                          * (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + l)
                              - *bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + l)))
                      : (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l])
                          + (1. / twoQ) * (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + l)
                              - *bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS
                                  + D3Q15::INVERSEDIRECTIONS[l])));
              }
            }
          }
        }

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_FINTERPOLATION_H */
