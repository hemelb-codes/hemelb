#include "lb/collisions/InletOutletWallCollision.h"
#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      InletOutletWallCollision::InletOutletWallCollision(distribn_t* iOutletDensityArray)
      {
        mBoundaryDensityArray = iOutletDensityArray;
      }

      void InletOutletWallCollision::Accept(Visitor* v,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl)
      {
        v->VisitInletOutletWall(this,
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
