
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H
#define HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * LBGK: This class implements the LBGK single-relaxation time kernel.
       */
      template<class KernelType, class LatticeType, class SecondLatticeType>
      class AdvectionDiffusionLBGK : public AdvectionDiffusionBaseKernel<KernelType, LatticeType, AdvectionDiffusionLBGK<KernelType, LatticeType, SecondLatticeType>, SecondLatticeType>
      {
        public:
          AdvectionDiffusionLBGK(AdvectionDiffusionInitParams& initParams)
          {
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<KernelType>& hydroVars, HydroVars<AdvectionDiffusionLBGK<KernelType, LatticeType, SecondLatticeType> >& advectionDiffusionVars, site_t index)
          {
            SecondLatticeType::CalculateDensityMomentumFEq(advectionDiffusionVars.f,
                                                           advectionDiffusionVars.density,
                                                           advectionDiffusionVars.momentum.x,
                                                           advectionDiffusionVars.momentum.y,
                                                           advectionDiffusionVars.momentum.z,
                                                           hydroVars.velocity.x,
                                                           hydroVars.velocity.y,
                                                           hydroVars.velocity.z,
                                                           advectionDiffusionVars.f_eq.f);

            for (unsigned int ii = 0; ii < SecondLatticeType::NUMVECTORS; ++ii)
            {
              advectionDiffusionVars.f_neq.f[ii] = advectionDiffusionVars.f[ii] - advectionDiffusionVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<KernelType>& hydroVars, HydroVars<AdvectionDiffusionLBGK<KernelType, LatticeType, SecondLatticeType> >& advectionDiffusionVars, site_t index)
          {
            SecondLatticeType::CalculateFeq(advectionDiffusionVars.density,
                                            hydroVars.velocity.x,
                                            hydroVars.velocity.y,
                                            hydroVars.velocity.z,
                                            advectionDiffusionVars.f_eq.f);

            for (unsigned int ii = 0; ii < SecondLatticeType::NUMVECTORS; ++ii)
            {
              advectionDiffusionVars.f_neq.f[ii] = advectionDiffusionVars.f[ii] - advectionDiffusionVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<AdvectionDiffusionLBGK<KernelType, LatticeType, SecondLatticeType> >& advectionDiffusionVars)
          {
            for (Direction direction = 0; direction < SecondLatticeType::NUMVECTORS; ++direction)
            {
              advectionDiffusionVars.SetFPostCollision(direction,
                                                       advectionDiffusionVars.f[direction]
                                                       + advectionDiffusionVars.f_neq.f[direction] * lbmParams->GetAdvectionDiffusionOmega());
            }
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H */
