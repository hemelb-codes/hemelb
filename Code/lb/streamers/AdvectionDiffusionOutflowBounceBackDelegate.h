
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWBOUNCEBACKDELEGATE_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWBOUNCEBACKDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/AdvectionDiffusionOutflowBaseStreamerDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class AdvectionDiffusionOutflowBounceBackDelegate : public AdvectionDiffusionOutflowBaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          AdvectionDiffusionOutflowBounceBackDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            collider(delegatorCollider)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const distribn_t &coupledV_x,
                                 const distribn_t &coupledV_y,
                                 const distribn_t &coupledV_z,
                                 const Direction& direction)
          {

            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

            distribn_t component = coupledV_x * LatticeType::CX[unstreamed] 
                                 + coupledV_y * LatticeType::CY[unstreamed]
                                 + coupledV_z * LatticeType::CZ[unstreamed];

            *latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
                = ((3.0 * component + 1) / (-3.0 * component + 1)) * hydroVars.GetFPostCollision()[direction];
          }
        protected:
          CollisionType& collider;
      };
    }
  }
}

#endif // HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWBOUNCEBACKDELEGATE_H
