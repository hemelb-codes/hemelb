#ifndef HEMELB_LB_COLLISIONS_IMPLGUOZHENGSHI_H
#define HEMELB_LB_COLLISIONS_IMPLGUOZHENGSHI_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      // Implementation of the Guo, Zheng, Shi boundary condition (2002).
      class ImplGuoZhengShi : public WallCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters* iLbmParams,
                            geometry::LatticeData* bLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters* iLbmParams,
                                    geometry::LatticeData* bLatDat,
                                    hemelb::vis::Control *iControl);

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLGUOZHENGSHI_H */
