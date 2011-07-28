#ifndef HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H
#define HEMELB_LB_COLLISIONS_STREAMANDCOLLIDE_H

#include "lb/collisions/CollisionVisitor.h"
#include "lb/collisions/implementations/Implementations.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      class StreamAndCollide : public CollisionVisitor
      {
        public:
              StreamAndCollide(tCollisionOperator* iCollisionOperator);

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
          tCollisionOperator* mCollisionOperator;

      }; /* End of StreamAndCollide definition */

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::StreamAndCollide(tCollisionOperator* iCollisionOperator)
      {
        mCollisionOperator = iCollisionOperator;
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::~StreamAndCollide()
      {

      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
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
                                                                   mCollisionOperator,
                                                                   iFirstIndex,
                                                                   iSiteCount,
                                                                   iLbmParams,
                                                                   bLatDat,
                                                                   iControl);
        }
        else
        {
          tInletOutletCollision::template DoStreamAndCollide<false>(mInletOutletCollision,
                                                                    mCollisionOperator,
                                                                    iFirstIndex,
                                                                    iSiteCount,
                                                                    iLbmParams,
                                                                    bLatDat,
                                                                    iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
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
                                                                       mCollisionOperator,
                                                                       iFirstIndex,
                                                                       iSiteCount,
                                                                       iLbmParams,
                                                                       bLatDat,
                                                                       iControl);
        }
        else
        {
          tInletOutletWallCollision::template DoStreamAndCollide<false>(mInletOutletWallCollision,
                                                                        mCollisionOperator,
                                                                        iFirstIndex,
                                                                        iSiteCount,
                                                                        iLbmParams,
                                                                        bLatDat,
                                                                        iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
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
                                                                mCollisionOperator,
                                                                iFirstIndex,
                                                                iSiteCount,
                                                                iLbmParams,
                                                                bLatDat,
                                                                iControl);
        }
        else
        {
          tMidFluidCollision::template DoStreamAndCollide<false>(mMidFluidCollision,
                                                                 mCollisionOperator,
                                                                 iFirstIndex,
                                                                 iSiteCount,
                                                                 iLbmParams,
                                                                 bLatDat,
                                                                 iControl);
        }
      }

      template<typename tMidFluidCollision, typename tWallCollision,
          typename tInletOutletCollision, typename tInletOutletWallCollision,
          typename tCollisionOperator>
      void StreamAndCollide<tMidFluidCollision, tWallCollision, tInletOutletCollision,
          tInletOutletWallCollision, tCollisionOperator>::VisitWall(WallCollision* mWallCollision,
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
                                                            mCollisionOperator,
                                                            iFirstIndex,
                                                            iSiteCount,
                                                            iLbmParams,
                                                            bLatDat,
                                                            iControl);
        }
        else
        {
          tWallCollision::template DoStreamAndCollide<false>(mWallCollision,
                                                             mCollisionOperator,
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
