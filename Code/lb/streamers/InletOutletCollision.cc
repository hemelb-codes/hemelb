#include "lb/streamers/InletOutletCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      InletOutletCollision::InletOutletCollision(boundaries::BoundaryValues* iBoundaryValues) :
        mBoundaryValues(iBoundaryValues)
      {

      }

      InletOutletCollision::~InletOutletCollision()
      {

      }

      distribn_t InletOutletCollision::getBoundaryDensityArray(const int index)
      {
        return mBoundaryValues->GetBoundaryDensity(index);
      }

      void InletOutletCollision::AcceptCollisionVisitor(collisions::CollisionVisitor* v,
                                                        const bool iDoRayTracing,
                                                        const site_t iFirstIndex,
                                                        const site_t iSiteCount,
                                                        const LbmParameters* iLbmParams,
                                                        geometry::LatticeData* bLatDat,
                                                        hemelb::vis::Control *iControl)
      {
        v->VisitInletOutlet(this,
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
