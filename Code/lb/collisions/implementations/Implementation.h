#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_IMPLEMENTATION_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_IMPLEMENTATION_H

#include "geometry/LatticeData.h"
#include "vis/Control.h"
#include "lb/LbmParameters.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        class Implementation
        {
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

        };
      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_IMPLEMENTATION_H */
