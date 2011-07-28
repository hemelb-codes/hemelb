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

      void CollisionVisitor::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                                  const bool iDoRayTracing,
                                                  const site_t iFirstIndex,
                                                  const site_t iSiteCount,
                                                  const LbmParameters* iLbmParams,
                                                  geometry::LatticeData* bLatDat,
                                                  hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                           const bool iDoRayTracing,
                                           const site_t iFirstIndex,
                                           const site_t iSiteCount,
                                           const LbmParameters* iLbmParams,
                                           geometry::LatticeData* bLatDat,
                                           hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void CollisionVisitor::VisitWall(WallCollision* mWallCollision,
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
