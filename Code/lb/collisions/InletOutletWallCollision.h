#ifndef HEMELB_LB_COLLISIONS_INLETOUTLETWALLCOLLISION_H
#define HEMELB_LB_COLLISIONS_INLETOUTLETWALLCOLLISION_H

#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class InletOutletWallCollision : public Collision
      {

        public:
          InletOutletWallCollision(BoundaryComms* iBoundaryComms);

          virtual ~InletOutletWallCollision();

          distribn_t getBoundaryDensityArray(const int index);

          virtual void AcceptCollisionVisitor(CollisionVisitor* v,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);

        private:
          BoundaryComms* mBoundaryComms;

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_INLETOUTLETWALLCOLLISION_H */
