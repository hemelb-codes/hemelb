
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Collision operator that uses an externally imposed density, zero velocity and relaxes
       * fully to equilibrium.
       */
      template<typename KernelType>
      class ZeroVelocityEquilibriumFixedDensity : public BaseCollision<ZeroVelocityEquilibriumFixedDensity<KernelType>,
          KernelType>
      {
        public:
          typedef KernelType CKernel;

          ZeroVelocityEquilibriumFixedDensity(kernels::InitParams& initParams) :
              BaseCollision<ZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>(), kernel(initParams), boundaryObject(initParams.boundaryObject)
          {
          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              const geometry::Site<geometry::LatticeData>& site)
          {
            hydroVars.density = boundaryObject->GetBoundaryDensity(site.GetIoletId());

            hydroVars.momentum = util::Vector3D<distribn_t>::Zero();

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

#endif /* HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
