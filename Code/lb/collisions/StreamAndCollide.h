#ifndef HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H
#define HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H

#include "lb/collisions/CollisionVisitor.h"
#include "lb/streamers/Implementations.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      class StreamAndCollide : public CollisionVisitor
      {
        public:
          virtual ~StreamAndCollide();

          virtual void VisitInletOutlet(streamers::InletOutletCollision* mInletOutletCollision,
                                        const bool iDoRayTracing,
                                        const site_t iFirstIndex,
                                        const site_t iSiteCount,
                                        const LbmParameters* iLbmParams,
                                        geometry::LatticeData* bLatDat,
                                        hemelb::vis::Control *iControl);

          virtual void VisitInletOutletWall(streamers::InletOutletWallCollision* mInletOutletWallCollision,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl);

          virtual void VisitMidFluid(streamers::MidFluidCollision* mMidFluidCollision,
                                     const bool iDoRayTracing,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl);

          virtual void VisitWall(streamers::WallCollision* mWallCollision,
                                 const bool iDoRayTracing,
                                 const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl);

      };
      /* End of StreamAndCollide definition */

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::~StreamAndCollide()
      {

      }

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitInletOutlet(streamers::InletOutletCollision* mInletOutletCollision,
                                                       const bool iDoRayTracing,
                                                       const site_t iFirstIndex,
                                                       const site_t iSiteCount,
                                                       const LbmParameters* iLbmParams,
                                                       geometry::LatticeData* bLatDat,
                                                       hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tInletOutletCollision::template DoStreamAndCollide<true>(mInletOutletCollision,
                                                                   iFirstIndex,
                                                                   iSiteCount,
                                                                   iLbmParams,
                                                                   bLatDat,
                                                                   iControl);
        }
        else
        {
          tInletOutletCollision::template DoStreamAndCollide<false>(mInletOutletCollision,
                                                                    iFirstIndex,
                                                                    iSiteCount,
                                                                    iLbmParams,
                                                                    bLatDat,
                                                                    iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitInletOutletWall(streamers::InletOutletWallCollision* mInletOutletWallCollision,
                                                           const bool iDoRayTracing,
                                                           const site_t iFirstIndex,
                                                           const site_t iSiteCount,
                                                           const LbmParameters* iLbmParams,
                                                           geometry::LatticeData* bLatDat,
                                                           hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tInletOutletWallCollision::template DoStreamAndCollide<true>(mInletOutletWallCollision,
                                                                       iFirstIndex,
                                                                       iSiteCount,
                                                                       iLbmParams,
                                                                       bLatDat,
                                                                       iControl);
        }
        else
        {
          tInletOutletWallCollision::template DoStreamAndCollide<false>(mInletOutletWallCollision,
                                                                        iFirstIndex,
                                                                        iSiteCount,
                                                                        iLbmParams,
                                                                        bLatDat,
                                                                        iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitMidFluid(streamers::MidFluidCollision* mMidFluidCollision,
                                                    const bool iDoRayTracing,
                                                    const site_t iFirstIndex,
                                                    const site_t iSiteCount,
                                                    const LbmParameters* iLbmParams,
                                                    geometry::LatticeData* bLatDat,
                                                    hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tMidFluidCollision::template DoStreamAndCollide<true>(mMidFluidCollision,
                                                                iFirstIndex,
                                                                iSiteCount,
                                                                iLbmParams,
                                                                bLatDat,
                                                                iControl);
        }
        else
        {
          tMidFluidCollision::template DoStreamAndCollide<false>(mMidFluidCollision,
                                                                 iFirstIndex,
                                                                 iSiteCount,
                                                                 iLbmParams,
                                                                 bLatDat,
                                                                 iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision, typename tInletOutletCollision,
          typename tInletOutletWallCollision>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision>::VisitWall(streamers::WallCollision* mWallCollision,
                                                const bool iDoRayTracing,
                                                const site_t iFirstIndex,
                                                const site_t iSiteCount,
                                                const LbmParameters* iLbmParams,
                                                geometry::LatticeData* bLatDat,
                                                hemelb::vis::Control *iControl)
      {
        if (iDoRayTracing)
        {
          tWallCollision::template DoStreamAndCollide<true>(mWallCollision,
                                                            iFirstIndex,
                                                            iSiteCount,
                                                            iLbmParams,
                                                            bLatDat,
                                                            iControl);
        }
        else
        {
          tWallCollision::template DoStreamAndCollide<false>(mWallCollision,
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

#endif /* HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H */
