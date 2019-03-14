
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONCOATINGCONCENTRATIONDELEGATE_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONCOATINGCONCENTRATIONDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class AdvectionDiffusionCoatingConcentrationDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          AdvectionDiffusionCoatingConcentrationDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            collider(delegatorCollider), bValues(initParams.advectionDiffusionBoundaryObject)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& direction)
          {
            int boundaryId = site.GetStentId() - 2;

            // Set the density at the "ghost" site to be the density of the stent.
            stents::StentCoating* stent =
                dynamic_cast<stents::StentCoating*>(bValues->GetLocalStent(boundaryId));
 
            distribn_t Dc = stent->GetCoatingDiffusivity();
            distribn_t lc = stent->GetCoatingThickness();
            distribn_t cs = stent->GetDensity(bValues->GetTimeStep());

            distribn_t k = (PI*PI*Dc) / (lc*lc);

            int N = 101;

            distribn_t c = 0.0;

            for(int i = 0; i < N; i++)
            {
             c = c + 2.0*((1.0)/((i+0.5)*(i+0.5)*PI*PI))*cs*std::exp(-(i+0.5)*(i+0.5)*k*bValues->GetTimeStep());
            }

            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

            *latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
                = -hydroVars.GetFPostCollision()[direction] + 2 * LatticeType::EQMWEIGHTS[direction] * c;
          }
        protected:
          CollisionType& collider;
          stents::BoundaryValues* bValues;
      };
    }
  }
}

#endif // HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONCOATINGCONCENTRATIONDELEGATE_H
