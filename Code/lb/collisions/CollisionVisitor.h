#ifndef HEMELB_LB_COLLISIONS_COLLISION_VISITOR_H
#define HEMELB_LB_COLLISIONS_COLLISION_VISITOR_H

#include "lb/collisions/Collisions.h"
#include "vis/Control.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"

#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      /*
       * The CollisionVisitor class separates the implementation of the streaming, collision and post step
       * from the object structure. This was done in order to make extensions to different collision
       * models easier via templating(impossible in the previous version as you cannot
       * have virtual template functions).
       */
      class CollisionVisitor
      {
        public:
          virtual ~CollisionVisitor();

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

        protected:
          template<bool tDoRayTracing>
          static void UpdateMinsAndMaxes(distribn_t iVx,
                                         distribn_t iVy,
                                         distribn_t iVz,
                                         const site_t iSiteIndex,
                                         const distribn_t* f_neq,
                                         const distribn_t iDensity,
                                         const geometry::LatticeData* iLatDat,
                                         const LbmParameters* iLbmParams,
                                         hemelb::vis::Control *iControl);
      };

      template<bool tDoRayTracing>
      void CollisionVisitor::UpdateMinsAndMaxes(distribn_t iVx,
                                                distribn_t iVy,
                                                distribn_t iVz,
                                                const site_t iSiteIndex,
                                                const distribn_t* f_neq,
                                                const distribn_t iDensity,
                                                const geometry::LatticeData* iLatDat,
                                                const LbmParameters* iLbmParams,
                                                hemelb::vis::Control *iControl)
      {
        if (tDoRayTracing)
        {
          distribn_t rtStress;

          if (iLbmParams->StressType == ShearStress)
          {
            if (iLatDat->GetNormalToWall(iSiteIndex)[0] > NO_VALUE)
            {
              rtStress = NO_VALUE;
            }
            else
            {
              D3Q15::CalculateShearStress(iDensity,
                                          f_neq,
                                          iLatDat->GetNormalToWall(iSiteIndex),
                                          rtStress,
                                          iLbmParams->StressParameter);
            }
          }
          else
          {
            D3Q15::CalculateVonMisesStress(f_neq, rtStress, iLbmParams->StressParameter);
          }

          // TODO: It'd be nice if the /iDensity were unnecessary.
          distribn_t lVelocity = sqrt(iVx * iVx + iVy * iVy + iVz * iVz) / iDensity;
          iControl->RegisterSite(iSiteIndex, iDensity, lVelocity, rtStress);
        }
      }

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_COLLISION_VISITOR_H */
