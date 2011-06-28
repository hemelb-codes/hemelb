#include "lb/collisions/MidFluidCollision.h"
#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      void MidFluidCollision::Accept(Visitor* v,
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
