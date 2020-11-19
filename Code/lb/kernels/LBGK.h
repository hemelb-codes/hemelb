// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_LBGK_H
#define HEMELB_LB_KERNELS_LBGK_H

#include <cstdlib>
#include <cmath>
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
      class LBGK : public BaseKernel<LBGK<LatticeType>, LatticeType>
      {
        public:
          LBGK(InitParams& initParams)
          {
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<LBGK<LatticeType> >& hydroVars,
                                                    site_t index)
          {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum.x,
                                                     hydroVars.momentum.y,
                                                     hydroVars.momentum.z,
                                                     hydroVars.velocity.x,
                                                     hydroVars.velocity.y,
                                                     hydroVars.velocity.z,
                                                     hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCalculateFeq(HydroVars<LBGK>& hydroVars, site_t index)
          {
            LatticeType::CalculateFeq(hydroVars.density,
                                      hydroVars.momentum.x,
                                      hydroVars.momentum.y,
                                      hydroVars.momentum.z,
                                      hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
              hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
          }

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<LBGK>& hydroVars)
          {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
              hydroVars.SetFPostCollision(direction,
                                          hydroVars.f[direction]
                                              + hydroVars.f_neq.f[direction]
                                                  * lbmParams->GetOmega());
          }

      };

    }
  }
}

#endif /* HEMELB_LB_KERNELS_LBGK_H */
