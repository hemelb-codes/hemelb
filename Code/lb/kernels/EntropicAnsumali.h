// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ENTROPICANSUMALI_H
#define HEMELB_LB_KERNELS_ENTROPICANSUMALI_H

#include "lb/concepts.h"
#include "lb/kernels/Entropic.h"
#include "lb/HydroVars.h"

namespace hemelb::lb
{
    // We have to declare this up here in order for it to be used as a template parameter in the
    // following declaration. Moving the template specialisation to the bottom of the file would
    // prevent it from being used as the HydroVars for this kernel.
    template<lattice_type LatticeType> class EntropicAnsumali;

    template<lattice_type LatticeType>
    struct HydroVars<EntropicAnsumali<LatticeType> > : public HydroVarsBase<LatticeType>
    {
        using HydroVarsBase<LatticeType>::HydroVarsBase;

        site_t index;
    };

    /**
     * EntropicBase: This class implements the entropic kernel as according to Anusmali et al, a modification to the standard
     * LBGK kernel which ensures the increase of entropy.
     */
    template<lattice_type L>
    class EntropicAnsumali : public EntropicBase<L>
    {
    public:
        using Base = EntropicBase<L>;
        using LatticeType = L;
        using VarsType = HydroVars<EntropicAnsumali>;

        /**
         * Constructor, passes parameters onto the base class.
         * @param initParams
         */
        EntropicAnsumali(InitParams& initParams) :
                Base(&initParams)
        {
        }

        /**
         * Calculates the density and momentum for the given f. Then calculates the
         * equilibrium distribution as described by Ansumali.
         * @param hydroVars
         * @param index The current lattice site index.
         */
        void CalculateDensityMomentumFeq(VarsType& hydroVars, site_t index)
        {
            hydroVars.index = index;
            LatticeType::CalculateDensityAndMomentum(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum);
            LatticeType::CalculateEntropicFeqAnsumali(hydroVars.density,
                                                      hydroVars.momentum,
                                                      hydroVars.f_eq);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }
        }

        /**
         * Calculates the equilibrium f distribution for the given density and momentum, as
         * described by Ansumali.
         * @param hydroVars
         * @param index The current lattice site index.
         */
        void CalculateFeq(VarsType& hydroVars, site_t index)
        {
            hydroVars.index = index;
            LatticeType::CalculateEntropicFeqAnsumali(hydroVars.density,
                                                      hydroVars.momentum,
                                                      hydroVars.f_eq);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }
        }
    };
}

#endif /* HEMELB_LB_KERNELS_ENTROPICANSUMALI_H */
