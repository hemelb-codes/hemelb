
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
      template<class LatticeType>
      class AdvectionDiffusionLBGK : public AdvectionDiffusionBaseKernel<AdvectionDiffusionLBGK<LatticeType>, LatticeType>
      {
        public:
          AdvectionDiffusionLBGK(InitParams& initParams)
          {
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<AdvectionDiffusionLBGK<LatticeType> >& hydroVars, lb::MacroscopicPropertyCache& propertyCache, site_t index)
          {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum.x,
                                                     hydroVars.momentum.y,
                                                     hydroVars.momentum.z,
                                                     propertyCache.velocityCache.Get(index).x,
                                                     propertyCache.velocityCache.Get(index).y,
                                                     propertyCache.velocityCache.Get(index).z,
                                                     hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<AdvectionDiffusionLBGK>& hydroVars, lb::MacroscopicPropertyCache& propertyCache, site_t index)
          {
            LatticeType::CalculateFeq(hydroVars.density,
                                      propertyCache.velocityCache.Get(index).x,
                                      propertyCache.velocityCache.Get(index).y,
                                      propertyCache.velocityCache.Get(index).z,
                                      hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<AdvectionDiffusionLBGK>& hydroVars)
          {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              hydroVars.SetFPostCollision(direction,
                                          hydroVars.f[direction]
                                              + hydroVars.f_neq.f[direction] * lbmParams->GetAdvectionDiffusionOmega());
            }
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_ADVECTIONDIFFUSIONLBGK_H */
