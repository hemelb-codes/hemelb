#include "lb/collisions/InletOutletWallCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      InletOutletWallCollision::InletOutletWallCollision(BoundaryComms* iBoundaryComms,
                                                         unsigned int iIOtype) :
        IOtype(iIOtype), mBoundaryComms(iBoundaryComms)
      {
        if (IOtype == INLET)
          mBoundaryDensityArray = mBoundaryComms->inlet_density;
        else if (IOtype == OUTLET)
          mBoundaryDensityArray = mBoundaryComms->outlet_density;
      }

      InletOutletWallCollision::~InletOutletWallCollision()
      {

      }

      distribn_t InletOutletWallCollision::getBoundaryDensityArray(const int index)
      {
        mBoundaryComms->WaitForComms(index, IOtype);

        return mBoundaryDensityArray[index];
      }

      void InletOutletWallCollision::AcceptCollisionVisitor(CollisionVisitor* v,
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
