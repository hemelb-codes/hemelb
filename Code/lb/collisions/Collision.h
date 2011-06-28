#ifndef COLLISION_H
#define COLLISION_H

#include "vis/Control.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
//#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class Visitor; // To deal with circular dependency of collision and visitor

      class Collision
      {
        public:
          virtual ~Collision();

          virtual void Accept(Visitor* v,
                              const bool iDoRayTracing,
                              const site_t iFirstIndex,
                              const site_t iSiteCount,
                              const LbmParameters* iLbmParams,
                              geometry::LatticeData* bLatDat,
                              hemelb::vis::Control *iControl) = 0;

      };

    }
  }
}

#endif /* COLLISION_H */
