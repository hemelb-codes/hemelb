#include "lb/collisions/WallCollision.h"
#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void WallCollision::Accept(Visitor* v,
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
