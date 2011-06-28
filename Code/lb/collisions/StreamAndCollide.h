#ifndef STREAMANDCOLLIDE_H
#define STREAMANDCOLLIDE_H

#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class StreamAndCollide : public Visitor
      {
        public:
          virtual ~StreamAndCollide();

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
          void FInterpolation(WallCollision* mWallCollision,
                              const site_t iFirstIndex,
                              const site_t iSiteCount,
                              const LbmParameters* iLbmParams,
                              geometry::LatticeData* bLatDat,
                              hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void GuoZhengShi(WallCollision* mWallCollision,
                           const site_t iFirstIndex,
                           const site_t iSiteCount,
                           const LbmParameters* iLbmParams,
                           geometry::LatticeData* bLatDat,
                           hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void NonZeroVelocityBoundaryDensity(InletOutletCollision* mInletOutletCollision,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void Regularised(WallCollision* mWallCollision,
                           const site_t iFirstIndex,
                           const site_t iSiteCount,
                           const LbmParameters* iLbmParams,
                           geometry::LatticeData* bLatDat,
                           hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void SimpleBounceBack(WallCollision* mWallCollision,
                                const site_t iFirstIndex,
                                const site_t iSiteCount,
                                const LbmParameters* iLbmParams,
                                geometry::LatticeData* bLatDat,
                                hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void SimpleCollideAndStream(MidFluidCollision* mMidFluidCollision,
                                      const site_t iFirstIndex,
                                      const site_t iSiteCount,
                                      const LbmParameters* iLbmParams,
                                      geometry::LatticeData* bLatDat,
                                      hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void ZeroVelocityBoundaryDensity(InletOutletWallCollision* mInletOutletWallCollision,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void ZeroVelocityEquilibrium(WallCollision* mWallCollision,
                                       const site_t iFirstIndex,
                                       const site_t iSiteCount,
                                       const LbmParameters* iLbmParams,
                                       geometry::LatticeData* bLatDat,
                                       hemelb::vis::Control *iControl);
      };

    }
  }
}

#endif /* STREAMANDCOLLIDE_H */
