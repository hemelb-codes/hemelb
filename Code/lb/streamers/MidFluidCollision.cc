#include "lb/streamers/MidFluidCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      MidFluidCollision::~MidFluidCollision()
      {

      }

      void MidFluidCollision::AcceptCollisionVisitor(collisions::CollisionVisitor* v,
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
