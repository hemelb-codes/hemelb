// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      template<typename KernelType>
      class NonZeroVelocityEquilibriumFixedDensity : public BaseCollision<
          NonZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          NonZeroVelocityEquilibriumFixedDensity(kernels::InitParams& initParams) :
              kernel(initParams), boundaryObject(initParams.boundaryObject)
          {

          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              const geometry::Site<geometry::LatticeData>& site)
          {
            CKernel::LatticeType::CalculateDensityAndMomentum(hydroVars.f,
                                                              hydroVars.density,
                                                              hydroVars.momentum.x,
                                                              hydroVars.momentum.y,
                                                              hydroVars.momentum.z);

            // Externally impose a density. Keep a record of the old one so we can scale the
            // momentum vector.
            distribn_t previousDensity = hydroVars.density;
            hydroVars.density = boundaryObject->GetBoundaryDensity(site.GetIoletId());

            hydroVars.momentum *= (hydroVars.density / previousDensity);

            kernel.CalculateFeq(hydroVars, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams, kernels::HydroVars<KernelType>& iHydroVars)
          {
            for (Direction direction = 0; direction < CKernel::LatticeType::NUMVECTORS; ++direction)
            {
              iHydroVars.SetFPostCollision(direction, iHydroVars.GetFEq()[direction]);
            }
          }

        private:
          KernelType kernel;
          iolets::BoundaryValues* boundaryObject;
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
