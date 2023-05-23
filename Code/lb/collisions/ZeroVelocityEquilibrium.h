// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUM_H

#include "lb/concepts.h"

namespace hemelb::lb
{
    /**
     * Collision operator that maintains the density, set velocity to 0 and fully relaxes
     * to equilibrium.
     */
    template<kernel_type K>
    class ZeroVelocityEquilibrium
    {
    public:
        using KernelType = K;
        using LatticeType = typename KernelType::LatticeType;
        using VarsType = typename KernelType::VarsType;

        ZeroVelocityEquilibrium(InitParams& initParams) : kernel(initParams)
        {
        }

        void CalculatePreCollision(VarsType& hydroVars,
                                   const geometry::Site<geometry::Domain>& site)
        {
            hydroVars.density = 0.0;

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.density += hydroVars.f[ii];
            }

            hydroVars.momentum = util::Vector3D<distribn_t>::Zero();

            kernel.CalculateFeq(hydroVars, site.GetIndex());
        }

        void Collide(const LbmParameters* lbmParams,
                     VarsType& iHydroVars)
        {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
                iHydroVars.SetFPostCollision(direction, iHydroVars.GetFEq()[direction]);
            }
        }

    private:
        KernelType kernel;
      };
}

#endif
