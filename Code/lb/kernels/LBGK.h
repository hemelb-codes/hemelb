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
#include "lb/LbmParameters.h"

namespace hemelb::lb::kernels
{
    /**
     * LBGK: This class implements the LBGK single-relaxation time kernel.
     */
    template<lattice_type LatticeType>
    class LBGK : public BaseKernel<LBGK<LatticeType>, LatticeType>
    {
    public:
        using Base = BaseKernel<LBGK<LatticeType>, LatticeType>;
        using typename Base::KHydroVars;
        using typename Base::LatticeType;

        LBGK(InitParams& initParams)
        {
        }

        void DoCalculateDensityMomentumFeq(KHydroVars& hydroVars,
                                           site_t index)
        {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum,
                                                     hydroVars.velocity,
                                                     hydroVars.f_eq);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }
        }

        void DoCalculateFeq(KHydroVars& hydroVars, site_t index)
        {
            LatticeType::CalculateFeq(hydroVars.density,
                                      hydroVars.momentum,
                                      hydroVars.f_eq);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }
        }

        void DoCollide(const LbmParameters* const lbmParams, KHydroVars& hydroVars)
        {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
                hydroVars.SetFPostCollision(direction,
                                            hydroVars.f[direction]
                                            + hydroVars.f_neq[direction]
                                              * lbmParams->GetOmega());
        }
    };
}
#endif /* HEMELB_LB_KERNELS_LBGK_H */
