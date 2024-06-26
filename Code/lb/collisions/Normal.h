// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_NORMAL_H
#define HEMELB_LB_COLLISIONS_NORMAL_H

#include "lb/concepts.h"

namespace hemelb::lb
{
    /**
     * Normal collisions - use a formula to relax the distribution towards equilibrium.
     */
    template<kernel_type K>
    class Normal
    {
    public:
        using KernelType = K;
        using LatticeType = typename KernelType::LatticeType;
        using VarsType = typename KernelType::VarsType;

        Normal(InitParams& initParams) :
                kernel(initParams)
        {
        }

        void CalculatePreCollision(VarsType& hydroVars,
                                   const geometry::Site<geometry::Domain>& site)
        {
            kernel.CalculateDensityMomentumFeq(hydroVars, site.GetIndex());
        }

        void Collide(const LbmParameters* lbmParams, VarsType& iHydroVars)
        {
            kernel.Collide(lbmParams, iHydroVars);
        }

        KernelType kernel;
    };

}

#endif /* HEMELB_LB_COLLISIONS_NORMAL_H */
