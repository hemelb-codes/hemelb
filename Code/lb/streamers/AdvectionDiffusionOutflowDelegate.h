
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWDELEGATE_H
#define HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/AdvectionDiffusionOutflowBaseStreamerDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      template<typename CollisionImpl>
      class AdvectionDiffusionOutflowDelegate : public AdvectionDiffusionOutflowBaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          AdvectionDiffusionOutflowDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            collider(delegatorCollider), iolet(*initParams.boundaryObject)
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
            int boundaryId = site.GetIoletId();

            kernels::HydroVars<typename CollisionType::CKernel> oldHydrovars(site.GetFOld<LatticeType> ());

            LatticeType::CalculateDensityAndMomentum(site.GetFOld<LatticeType> (), oldHydrovars.density, oldHydrovars.momentum.x, oldHydrovars.momentum.y, oldHydrovars.momentum.z);

            util::Vector3D<float> ioletNormal = iolet.GetLocalIolet(boundaryId)->GetNormal();

            distribn_t component = coupledV_x * ioletNormal.x 
                                 + coupledV_y * ioletNormal.y
                                 + coupledV_z * ioletNormal.z;

            component *= oldHydrovars.density;

            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

            *latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
                = hydroVars.GetFPostCollision()[direction] + component;
          }
        protected:
          CollisionType& collider;
          iolets::BoundaryValues& iolet;
      };
    }
  }
}

#endif // HEMELB_LB_STREAMERS_ADVECTIONDIFFUSIONOUTFLOWDELEGATE_H
