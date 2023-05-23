// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/concepts.h"

namespace hemelb::lb
{
    /**
     * Collision operator that uses an externally imposed density, zero velocity and relaxes
     * fully to equilibrium.
     */
    template<kernel_type K>
    class ZeroVelocityEquilibriumFixedDensity
    {
    public:
        using KernelType = K;
        using LatticeType = typename KernelType::LatticeType;
        using VarsType = typename KernelType::VarsType;

        ZeroVelocityEquilibriumFixedDensity(InitParams& initParams) :
                kernel(initParams), boundaryObject(initParams.boundaryObject)
        {
        }

        void CalculatePreCollision(VarsType& hydroVars,
                                   const geometry::Site<geometry::Domain>& site)
        {
            hydroVars.density = boundaryObject->GetBoundaryDensity(site.GetIoletId());

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
        BoundaryValues* boundaryObject;

    };
}

#endif
