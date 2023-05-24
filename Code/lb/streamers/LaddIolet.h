// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_LADDIOLET_H
#define HEMELB_LB_STREAMERS_LADDIOLET_H

#include "lb/concepts.h"
#include "lb/streamers/SimpleBounceBack.h"

namespace hemelb::lb
{

    template<collision_type C>
    class LaddIoletLink : public BounceBackLink<C>
    {
    public:
        using CollisionType = C;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;

        LaddIoletLink(CollisionType& delegatorCollider, InitParams& initParams) :
                BounceBackLink<CollisionType>(delegatorCollider, initParams),
                bValues(initParams.boundaryObject)
        {
        }

        void StreamLink(const LbmParameters* lbmParams,
                               geometry::FieldData& latticeData,
                               const geometry::Site<geometry::Domain>& site,
                               VarsType& hydroVars,
                               const Direction& ii)
        {
            // Translating from Ladd, J. Fluid Mech. "Numerical simulations
            // of particulate suspensions via a discretized Boltzmann
            // equation. Part 1. Theoretical foundation", 1994
            // Eq (3.2) -- simple bounce-back -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*)
            // Eq (3.3) --- modified BB -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*) - 2 a1_i \rho u . c_i
            // where u is the velocity of the boundary half way along the
            // link and a1_i = w_1 / cs2

            int boundaryId = site.GetIoletId();
            auto iolet = dynamic_cast<InOutLetVelocity*>(bValues->GetLocalIolet(boundaryId));
            LatticePosition sitePos(site.GetGlobalSiteCoords());

            LatticePosition halfWay(sitePos);
            halfWay += 0.5 * LatticeType::VECTORS[ii];

            LatticeVelocity wallMom(iolet->GetVelocity(halfWay, bValues->GetTimeStep()));
            //TODO: Add site.GetGlobalSiteCoords() as a first argument?

            if (LatticeType::IsLatticeCompressible())
            {
                wallMom *= hydroVars.density;
            }

            distribn_t correction = 2. * LatticeType::EQMWEIGHTS[ii]
                                    * Dot(wallMom, LatticeType::VECTORS[ii]) / Cs2;

            * (latticeData.GetFNew(BounceBackLink<CollisionType>::GetBBIndex(site.GetIndex(),
                                                                             ii))) =
                    hydroVars.GetFPostCollision()[ii] - correction;
        }
    private:
        BoundaryValues* bValues;
    };

}

#endif
