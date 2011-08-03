#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      CollisionVisitor::~CollisionVisitor()
      {

      }

      void CollisionVisitor::VisitInletOutlet(streamers::InletOutletCollision* mInletOutletCollision,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitInletOutletWall(streamers::InletOutletWallCollision* mInletOutletWallCollision,
                                                  const bool iDoRayTracing,
                                                  const site_t iFirstIndex,
                                                  const site_t iSiteCount,
                                                  const LbmParameters* iLbmParams,
                                                  geometry::LatticeData* bLatDat,
                                                  hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitMidFluid(streamers::MidFluidCollision* mMidFluidCollision,
                                           const bool iDoRayTracing,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitWall(streamers::WallCollision* mWallCollision,
                                       const bool iDoRayTracing,
                                       const site_t iFirstIndex,
                                       const site_t iSiteCount,
                                       const LbmParameters* iLbmParams,
                                       geometry::LatticeData* bLatDat,
                                       hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

    }
  }
}
