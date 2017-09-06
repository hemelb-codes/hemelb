
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONVESSELWALLABSORPTIONDELEGATE_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONVESSELWALLABSORPTIONDELEGATE_H

#include "lb/streamers/BaseStreamerDelegate.h"
#include "lb/streamers/SimpleCollideAndStream.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class AdvectionDiffusionVesselWallAbsorptionDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          static inline site_t GetBBIndex(site_t siteIndex, int direction)
          {
            return (siteIndex * LatticeType::NUMVECTORS) + LatticeType::INVERSEDIRECTIONS[direction];
          }

          AdvectionDiffusionVesselWallAbsorptionDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            // Propagate the outgoing post-collisional f into the opposite direction.
            * (latticeData->GetFNew(GetBBIndex(site.GetIndex(), direction))) = ((-3 * lbmParams->GetAdvectionDiffusionAlpha() + 1) / (3 * lbmParams->GetAdvectionDiffusionAlpha() + 1))
                                                                                * hydroVars.GetFPostCollision()[direction];
          }

      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONVESSELWALLABSORPTIONDELEGATE_H */
