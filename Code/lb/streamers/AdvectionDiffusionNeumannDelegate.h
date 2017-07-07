
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONNEUMANNDELEGATE_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONNEUMANNDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"
#include "lb/stents/StentFlux.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class AdvectionDiffusionNeumannDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          AdvectionDiffusionNeumannDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            collider(delegatorCollider), bValues(initParams.boundaryObject)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            int boundaryId = site.GetIoletId();

            // Set the density at the "ghost" site to be the density of the stent.
            stents::StentFlux* stent =
                dynamic_cast<stents::StentFlux*>(bValues->GetLocalStent(boundaryId));
 
            LatticeSpeed stentwallFlux(stent->GetFlux(bValues->GetTimeStep()));

            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

            *latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
                = hydroVars.GetFPostCollision()[direction]+stentwallFlux;
          }
        protected:
          CollisionType& collider;
          stents::BoundaryValues* bValues;
      };
    }
  }
}

#endif // HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONNEUMANNDELEGATE_H
