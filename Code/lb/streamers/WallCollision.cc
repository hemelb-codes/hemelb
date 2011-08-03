#include "lb/streamers/WallCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      void WallCollision::AcceptCollisionVisitor(collisions::CollisionVisitor* v,
                                                 const bool iDoRayTracing,
                                                 const site_t iFirstIndex,
                                                 const site_t iSiteCount,
                                                 const LbmParameters* iLbmParams,
                                                 geometry::LatticeData* bLatDat,
                                                 hemelb::vis::Control *iControl)
      {
        v->VisitWall(this, iDoRayTracing, iFirstIndex, iSiteCount, iLbmParams, bLatDat, iControl);
      }

    }
  }
}
