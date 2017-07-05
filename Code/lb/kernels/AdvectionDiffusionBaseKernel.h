
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONBASEKERNEL_H
#define HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONBASEKERNEL_H

#include <cstdlib>
#include "constants.h"
#include "lb/kernels/BaseKernel.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {

      template<typename KernelImpl, typename LatticeImpl>
      class AdvectionDiffusionBaseKernel
      {
        public:
          typedef HydroVars<KernelImpl> KHydroVars;
          typedef LatticeImpl LatticeType;

          inline void CalculateDensityMomentumFeq(KHydroVars& hydroVars, KHydroVars& advectionDiffusionVars, site_t index)
          {
            static_cast<KernelImpl*> (this)->DoCalculateDensityMomentumFeq(hydroVars, advectionDiffusionVars, index);
          }

          inline void CalculateFeq(KHydroVars& hydroVars, KHydroVars& advectionDiffusionVars, site_t index)
          {
            static_cast<KernelImpl*> (this)->DoCalculateFeq(hydroVars, advectionDiffusionVars, index);
          }

          inline void Collide(const LbmParameters* lbmParams, KHydroVars& hydroVars)
          {
            static_cast<KernelImpl*> (this)->DoCollide(lbmParams, hydroVars);
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONBASEKERNEL_H */
