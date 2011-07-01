#ifndef HEMELB_LB_COLLISIONS_POSTSTEP_H
#define HEMELB_LB_COLLISIONS_POSTSTEP_H

#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      template<bool tDoEntropic>
      class PostStep : public CollisionVisitor
      {
        public:
          virtual ~PostStep();

          virtual void VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                        const bool iDoRayTracing,
                                        const site_t iFirstIndex,
                                        const site_t iSiteCount,
                                        const LbmParameters* iLbmParams,
                                        geometry::LatticeData* bLatDat,
                                        hemelb::vis::Control *iControl);

          virtual void VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl);

          virtual void VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                     const bool iDoRayTracing,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl);

          virtual void VisitWall(WallCollision* mWallCollision,
                                 const bool iDoRayTracing,
                                 const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void FInterpolationPostStep(WallCollision* mWallCollision,
                                      const site_t iFirstIndex,
                                      const site_t iSiteCount,
                                      const LbmParameters* iLbmParams,
                                      geometry::LatticeData* bLatDat,
                                      hemelb::vis::Control *iControl);
      };

      template<bool tDoEntropic>
      PostStep<tDoEntropic>::~PostStep()
      {

      }

      template<bool tDoEntropic>
      void PostStep<tDoEntropic>::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                                   const bool iDoRayTracing,
                                                   const site_t iFirstIndex,
                                                   const site_t iSiteCount,
                                                   const LbmParameters* iLbmParams,
                                                   geometry::LatticeData* bLatDat,
                                                   hemelb::vis::Control *iControl)
      {

      }

      template<bool tDoEntropic>
      void PostStep<tDoEntropic>::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                                       const bool iDoRayTracing,
                                                       const site_t iFirstIndex,
                                                       const site_t iSiteCount,
                                                       const LbmParameters* iLbmParams,
                                                       geometry::LatticeData* bLatDat,
                                                       hemelb::vis::Control *iControl)
      {

      }

      template<bool tDoEntropic>
      void PostStep<tDoEntropic>::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                                const bool iDoRayTracing,
                                                const site_t iFirstIndex,
                                                const site_t iSiteCount,
                                                const LbmParameters* iLbmParams,
                                                geometry::LatticeData* bLatDat,
                                                hemelb::vis::Control *iControl)
      {

      }

      template<bool tDoEntropic>
      void PostStep<tDoEntropic>::VisitWall(WallCollision* mWallCollision,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl)
      {

      }

      template<bool tDoEntropic>
      template<bool tDoRayTracing>
      void PostStep<tDoEntropic>::FInterpolationPostStep(WallCollision* mWallCollision,
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

#endif /* HEMELB_LB_COLLISIONS_POSTSTEP_H */
