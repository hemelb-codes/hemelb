#include "lb/collisions/PostStep.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      PostStep::~PostStep()
      {

      }

      void PostStep::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                      const bool iDoRayTracing,
                                      const site_t iFirstIndex,
                                      const site_t iSiteCount,
                                      const LbmParameters* iLbmParams,
                                      geometry::LatticeData* bLatDat,
                                      hemelb::vis::Control *iControl)
      {

      }

      void PostStep::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                          const bool iDoRayTracing,
                                          const site_t iFirstIndex,
                                          const site_t iSiteCount,
                                          const LbmParameters* iLbmParams,
                                          geometry::LatticeData* bLatDat,
                                          hemelb::vis::Control *iControl)
      {

      }

      void PostStep::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                   const bool iDoRayTracing,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl)
      {

      }

      void PostStep::VisitWall(WallCollision* mWallCollision,
                               const bool iDoRayTracing,
                               const site_t iFirstIndex,
                               const site_t iSiteCount,
                               const LbmParameters* iLbmParams,
                               geometry::LatticeData* bLatDat,
                               hemelb::vis::Control *iControl)
      {

      }

      template<bool tDoRayTracing>
      void PostStep::FInterpolationPostStep(WallCollision* mWallCollision,
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
