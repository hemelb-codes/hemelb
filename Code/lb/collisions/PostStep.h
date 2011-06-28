#ifndef POSTSTEP_H
#define POSTSTEP_H

#include "lb/collisions/Visitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class PostStep : public Visitor
      {
        public:
          virtual ~PostStep();

          virtual void VisitInletOutlet(InletOutletCollision* mInletOutletCollision,
                                        const bool iDoRayTracing,
                                        const site_t iFirstIndex,
                                        const site_t iSiteCount,
                                        const LbmParameters* iLbmParams,
                                        geometry::LatticeData* bLatDat,
                                        hemelb::vis::Control *iControl);

          virtual void VisitInletOutletWall(InletOutletWallCollision* mInletOutletWallCollision,
                                            const bool iDoRayTracing,
                                            const site_t iFirstIndex,
                                            const site_t iSiteCount,
                                            const LbmParameters* iLbmParams,
                                            geometry::LatticeData* bLatDat,
                                            hemelb::vis::Control *iControl);

          virtual void VisitMidFluid(MidFluidCollision* mMidFluidCollision,
                                     const bool iDoRayTracing,
                                     const site_t iFirstIndex,
                                     const site_t iSiteCount,
                                     const LbmParameters* iLbmParams,
                                     geometry::LatticeData* bLatDat,
                                     hemelb::vis::Control *iControl);

          virtual void VisitWall(WallCollision* mWallCollision,
                                 const bool iDoRayTracing,
                                 const site_t iFirstIndex,
                                 const site_t iSiteCount,
                                 const LbmParameters* iLbmParams,
                                 geometry::LatticeData* bLatDat,
                                 hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void FInterpolationPostStep(WallCollision* mWallCollision,
                                      const site_t iFirstIndex,
                                      const site_t iSiteCount,
                                      const LbmParameters* iLbmParams,
                                      geometry::LatticeData* bLatDat,
                                      hemelb::vis::Control *iControl);
      };

    }
  }
}

#endif /* POSTSTEP_H */
