#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void Visitor::VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                     const bool iDoRayTracing,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Visitor::VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                         const bool iDoRayTracing,
                                         const site_t iFirstIndex,
                                         const site_t iSiteCount,
                                         const LbmParameters* iLbmParams,
                                         geometry::LatticeData* bLatDat,
                                         hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Visitor::VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                  const bool iDoRayTracing,
                                  const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Visitor::VisitWall(WallCollision* mWallCollision,
                              const bool iDoRayTracing,
                              const site_t iFirstIndex,
                              const site_t iSiteCount,
                              const LbmParameters* iLbmParams,
                              geometry::LatticeData* bLatDat,
                              hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      Visitor::Visitor()
      {

      }

    }
  }
}
