// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_ENTROPICCHIK_H
#define HEMELB_LB_KERNELS_ENTROPICCHIK_H

#include <cstdlib>

#include "lb/kernels/Entropic.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb::lb::kernels
{
    // We have to declare this up here in order for it to be used as a template parameter in the
    // following declaration. Moving the template specialisation to the bottom of the file would
    // prevent it from being used as the HydroVars for this kernel.
    template<lattice_type> class EntropicChik;

    template<lattice_type LatticeType>
    struct HydroVars<EntropicChik<LatticeType> > : public HydroVarsBase<LatticeType>
    {
        using HydroVarsBase<LatticeType>::HydroVarsBase;

        site_t index;
    };

    /**
     * EntropicChik: This class implements the entropic kernel, as per Chitakamarla et al.
     */
    template<lattice_type LatticeType>
    class EntropicChik : public BaseKernel<EntropicChik<LatticeType>, LatticeType>,
                         public Entropic<LatticeType>
    {
    public:
        /**
         * Constructor, passes parameters onto the base class.
         * @param initParams
         */
        EntropicChik(InitParams& initParams) :
                Entropic<LatticeType>(&initParams)
        {
        }

        /**
         * Calculates the density and momentum for the given f. Then calculates the
         * equilibrium distribution as described by Chikatamarla.
         * @param hydroVars
         * @param index, the current lattice site index.
         */
        inline void DoCalculateDensityMomentumFeq(
                HydroVars<EntropicChik<LatticeType> >& hydroVars, site_t index)
        {
            hydroVars.index = index;
            LatticeType::CalculateDensityAndMomentum(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum);
            LatticeType::CalculateEntropicFeqChik(hydroVars.density,
                                                  hydroVars.momentum,
                                                  hydroVars.f_eq);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq[ii] = hydroVars.f[ii] - hydroVars.f_eq[ii];
            }
        }

        /**
         * Calculates the equilibrium f distribution for the given density and momentum, as
         * described by Chikatamarla.
         * @param hydroVars
         * @param index The current lattice site index.
         */
        inline void DoCalculateFeq(HydroVars<EntropicChik<LatticeType> >& hydroVars, site_t index)
        {
            hydroVars.index = index;
            LatticeType::CalculateEntropicFeqChik(hydroVars.density,
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
