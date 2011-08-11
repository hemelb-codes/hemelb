#include "lb/collisions/InletOutletCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      InletOutletCollision::InletOutletCollision(BoundaryComms* iBoundaryComms,
                                                 unsigned int iIOtype) :
        IOtype(iIOtype), mBoundaryComms(iBoundaryComms)
      {
        if (IOtype == INLET)
          mBoundaryDensityArray = mBoundaryComms->inlet_density;
        else if (IOtype == OUTLET)
          mBoundaryDensityArray = mBoundaryComms->outlet_density;
      }

      InletOutletCollision::~InletOutletCollision()
      {

      }

      distribn_t InletOutletCollision::getBoundaryDensityArray(const int index)
      {
        mBoundaryComms->WaitForComms(index, IOtype);

        return mBoundaryDensityArray[index];
      }

      void InletOutletCollision::AcceptCollisionVisitor(CollisionVisitor* v,
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
