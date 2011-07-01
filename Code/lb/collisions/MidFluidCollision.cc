#include "lb/collisions/MidFluidCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void MidFluidCollision::AcceptCollisionVisitor(CollisionVisitor* v,
                                                     const bool iDoRayTracing,
                                                     const site_t iFirstIndex,
                                                     const site_t iSiteCount,
                                                     const LbmParameters* iLbmParams,
                                                     geometry::LatticeData* bLatDat,
                                                     hemelb::vis::Control *iControl)
      {
        v->VisitMidFluid(this,
                         iDoRayTracing,
                         iFirstIndex,
                         iSiteCount,
                         iLbmParams,
                         bLatDat,
                         iControl);
      }

    }
  }
}
