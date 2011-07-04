#ifndef HEMELB_LB_COLLISIONS_POSTSTEP_H
#define HEMELB_LB_COLLISIONS_POSTSTEP_H

#include "lb/collisions/CollisionVisitor.h"
#include "lb/collisions/implementations/Implementations.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
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

      };

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
      PostStep<tMidFluidCollision, tWallCollision, tInletOutletCollision, tInletOutletWallCollision>::~PostStep()
      {

      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
      void PostStep<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                                       const bool iDoRayTracing,
                                                       const site_t iFirstIndex,
                                                       const site_t iSiteCount,
                                                       const LbmParameters* iLbmParams,
                                                       geometry::LatticeData* bLatDat,
                                                       hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tInletOutletCollision::template DoPostStep<true>(mInletOutletCollision,
                                                           iFirstIndex,
                                                           iSiteCount,
                                                           iLbmParams,
                                                           bLatDat,
                                                           iControl);
        }
        else
        {
          tInletOutletCollision::template DoPostStep<false>(mInletOutletCollision,
                                                            iFirstIndex,
                                                            iSiteCount,
                                                            iLbmParams,
                                                            bLatDat,
                                                            iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
      void PostStep<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                                           const bool iDoRayTracing,
                                                           const site_t iFirstIndex,
                                                           const site_t iSiteCount,
                                                           const LbmParameters* iLbmParams,
                                                           geometry::LatticeData* bLatDat,
                                                           hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tInletOutletWallCollision::template DoPostStep<true>(mInletOutletWallCollision,
                                                               iFirstIndex,
                                                               iSiteCount,
                                                               iLbmParams,
                                                               bLatDat,
                                                               iControl);
        }
        else
        {
          tInletOutletWallCollision::template DoPostStep<false>(mInletOutletWallCollision,
                                                                iFirstIndex,
                                                                iSiteCount,
                                                                iLbmParams,
                                                                bLatDat,
                                                                iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
      void PostStep<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                                    const bool iDoRayTracing,
                                                    const site_t iFirstIndex,
                                                    const site_t iSiteCount,
                                                    const LbmParameters* iLbmParams,
                                                    geometry::LatticeData* bLatDat,
                                                    hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tMidFluidCollision::template DoPostStep<true>(mMidFluidCollision,
                                                        iFirstIndex,
                                                        iSiteCount,
                                                        iLbmParams,
                                                        bLatDat,
                                                        iControl);
        }
        else
        {
          tMidFluidCollision::template DoPostStep<false>(mMidFluidCollision,
                                                         iFirstIndex,
                                                         iSiteCount,
                                                         iLbmParams,
                                                         bLatDat,
                                                         iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision>
      void PostStep<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitWall(WallCollision* mWallCollision,
                                                const bool iDoRayTracing,
                                                const site_t iFirstIndex,
                                                const site_t iSiteCount,
                                                const LbmParameters* iLbmParams,
                                                geometry::LatticeData* bLatDat,
                                                hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tWallCollision::template DoPostStep<true>(mWallCollision,
                                                    iFirstIndex,
                                                    iSiteCount,
                                                    iLbmParams,
                                                    bLatDat,
                                                    iControl);
        }
        else
        {
          tWallCollision::template DoPostStep<false>(mWallCollision,
                                                     iFirstIndex,
                                                     iSiteCount,
                                                     iLbmParams,
                                                     bLatDat,
                                                     iControl);
        }
      }

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_POSTSTEP_H */
