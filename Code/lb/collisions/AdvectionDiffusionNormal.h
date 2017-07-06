
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_NORMAL_H
#define HEMELB_LB_COLLISIONS_NORMAL_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Normal collisions - use a formula to relax the distribution towards equilibrium.
       */
      template<typename KernelType, typename SecondKernelType>
      class AdvectionDiffusionNormal : public AdvectionDiffusionBaseCollision<AdvectionDiffusionNormal<KernelType, SecondKernelType>, KernelType, SecondKernelType>
      {
        public:
          typedef KernelType CKernel;
          typedef SecondKernelType SecondCKernel;

          Normal(kernels::AdvectionDiffusionInitParams& initParams) :
              kernel(initParams)
          {
          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              kernels::HydroVars<SecondKernelType>& advectionDiffusionVars,
                                              const geometry::Site<geometry::LatticeData>& site)
          {
            kernel.CalculateDensityMomentumFeq(hydroVars, advectionDiffusionVars, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams,
                                kernels::HydroVars<SecondKernelType>& iAdvectionDiffusionVars)
          {
            kernel.Collide(lbmParams, iAdvectionDiffusionVars);
          }


          SecondKernelType kernel;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_NORMAL_H */
