
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H
#define HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/AdvectionDiffusionBaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * LBGK: This class implements the LBGK single-relaxation time kernel.
       */
      template<class LatticeType, class AdvectionDiffusionLatticeType>
      class AdvectionDiffusionLBGK : public AdvectionDiffusionBaseKernel<AdvectionDiffusionLBGK<AdvectionDiffusionLatticeType>, AdvectionDiffusionLatticeType>
      {
        public:
          AdvectionDiffusionLBGK(InitParams& initParams)
          {
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<LBGK<LatticeType> >& hydroVars, HydroVars<AdvectionDiffusionLBGK<AdvectionDiffusionLatticeType> >& advectionDiffusionHydroVars, site_t index)
          {
            AdvectionDiffusionLatticeType::CalculateDensityMomentumFEq(advectionDiffusionHydroVars.f,
                                                                       advectionDiffusionHydroVars.density,
                                                                       advectionDiffusionHydroVars.momentum.x,
                                                                       advectionDiffusionHydroVars.momentum.y,
                                                                       advectionDiffusionHydroVars.momentum.z,
                                                                       hydroVars.velocity.x,
                                                                       hydroVars.velocity.y,
                                                                       hydroVars.velocity.z,
                                                                       advectionDiffusionHydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              advectionDiffusionHydroVars.f_neq.f[ii] = advectionDiffusionHydroVars.f[ii] - advectionDiffusionHydroVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<LBGK<LatticeType> >& hydroVars, HydroVars<AdvectionDiffusionLBGK>& advectionDiffusionHydroVars, site_t index)
          {
            AdvectionDiffusionLatticeType::CalculateFeq(advectionDiffusionHydroVars.density,
                                                        hydroVars.momentum.x/hydroVars.density,
                                                        hydroVars.momentum.y/hydroVars.density,
                                                        hydroVars.momentum.z/hydroVars.density,
                                                        advectionDiffusionHydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < AdvectionDiffusionLatticeType::NUMVECTORS; ++ii)
            {
              advectionDiffusionHydroVars.f_neq.f[ii] = advectionDiffusionHydroVars.f[ii] - advectionDiffusionHydroVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<AdvectionDiffusionLBGK>& advectionDiffusionHydroVars)
          {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              advectionDiffusionHydroVars.SetFPostCollision(direction,
                                                            advectionDiffusionHydroVars.f[direction]
                                                            + advectionDiffusionHydroVars.f_neq.f[direction] * lbmParams->GetAdvectionDiffusionOmega());
            }
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H */
