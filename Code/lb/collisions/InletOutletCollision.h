#ifndef HEMELB_LB_COLLISIONS_INLETOUTLETCOLLISION_H
#define HEMELB_LB_COLLISIONS_INLETOUTLETCOLLISION_H

#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class InletOutletCollision : public Collision
      {
        public:
          InletOutletCollision(BoundaryComms* iBoundaryComms, unsigned int iIOtype);

          ~InletOutletCollision();

          distribn_t getBoundaryDensityArray(const int index);

          virtual void AcceptCollisionVisitor(CollisionVisitor* v,
                                              const bool iDoRayTracing,
                                              const site_t iFirstIndex,
                                              const site_t iSiteCount,
                                              const LbmParameters* iLbmParams,
                                              geometry::LatticeData* bLatDat,
                                              hemelb::vis::Control *iControl);

          const unsigned int IOtype;

        private:
          distribn_t* mBoundaryDensityArray;
          BoundaryComms* mBoundaryComms;
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_INLETOUTLETCOLLISION_H */
