// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_COMMON_H
#define HEMELB_LB_STREAMERS_COMMON_H

#include <cmath>

#include "geometry/Domain.h"
#include "lb/concepts.h"
#include "lb/HydroVars.h"
#include "lb/LbmParameters.h"
#include "lb/MacroscopicPropertyCache.h"

namespace hemelb::lb
{
    /// Update the property cache.
    template<lattice_type LatticeType>
    void UpdateCachePostCollision(
            const geometry::Site<geometry::Domain>& site,
            const HydroVarsBase<LatticeType>& hydroVars, const LbmParameters* lbmParams,
            lb::MacroscopicPropertyCache& propertyCache)
    {
            if (propertyCache.densityCache.RequiresRefresh())
            {
              propertyCache.densityCache.Put(site.GetIndex(), hydroVars.density);
            }

            if (propertyCache.velocityCache.RequiresRefresh())
            {
              propertyCache.velocityCache.Put(site.GetIndex(), hydroVars.velocity);
            }

            if (propertyCache.wallShearStressMagnitudeCache.RequiresRefresh())
            {
              distribn_t stress;

              if (!site.IsWall())
              {
                stress = NO_VALUE;
              }
              else
              {
                LatticeType::CalculateWallShearStressMagnitude(hydroVars.density,
                                                               hydroVars.GetFNeq(),
                                                               site.GetWallNormal(),
                                                               stress,
                                                               lbmParams->GetStressParameter());
              }

              propertyCache.wallShearStressMagnitudeCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.vonMisesStressCache.RequiresRefresh())
            {
              distribn_t stress;
              LatticeType::CalculateVonMisesStress(hydroVars.GetFNeq(),
                                                                                         stress,
                                                                                         lbmParams->GetStressParameter());

              propertyCache.vonMisesStressCache.Put(site.GetIndex(), stress);
            }

            if (propertyCache.shearRateCache.RequiresRefresh())
            {
              distribn_t shear_rate =
                  LatticeType::CalculateShearRate(hydroVars.tau,
                                                                                        hydroVars.GetFNeq(),
                                                                                        hydroVars.density);

              propertyCache.shearRateCache.Put(site.GetIndex(), shear_rate);
            }

            if (propertyCache.stressTensorCache.RequiresRefresh())
            {
              util::Matrix3D stressTensor;
              LatticeType::CalculateStressTensor(hydroVars.density,
                                                                                       hydroVars.tau,
                                                                                       hydroVars.GetFNeq(),
                                                                                       stressTensor);

              propertyCache.stressTensorCache.Put(site.GetIndex(), stressTensor);

            }

            if (propertyCache.tractionCache.RequiresRefresh())
            {
              util::Vector3D<LatticeStress> tractionOnAPoint(0);

              /*
               * Wall normals are only available at the sites marked as being at the domain edge.
               * For the sites in the fluid bulk, the traction vector will be 0.
               */
              if (site.IsWall())
              {
                LatticeType::CalculateTractionOnAPoint(hydroVars.density,
                                                       hydroVars.tau,
                                                       hydroVars.GetFNeq(),
                                                       site.GetWallNormal(),
                                                       tractionOnAPoint);
              }

              propertyCache.tractionCache.Put(site.GetIndex(), tractionOnAPoint);

            }

            if (propertyCache.tangentialProjectionTractionCache.RequiresRefresh())
            {
              util::Vector3D<LatticeStress> tangentialProjectionTractionOnAPoint(0);

              /*
               * Wall normals are only available at the sites marked as being at the domain edge.
               * For the sites in the fluid bulk, the traction vector will be 0.
               */
              if (site.IsWall())
              {
                LatticeType::CalculateTangentialProjectionTraction(hydroVars.density,
                                                                   hydroVars.tau,
                                                                   hydroVars.GetFNeq(),
                                                                   site.GetWallNormal(),
                                                                   tangentialProjectionTractionOnAPoint);
              }

              propertyCache.tangentialProjectionTractionCache.Put(site.GetIndex(),
                                                                  tangentialProjectionTractionOnAPoint);

            }
          }

    /**
     * Null implementation of an iolet link delegate.
     */
    template<typename C>
    struct NullLink
    {
        using CollisionType = C;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;
        NullLink(CollisionType& collider, InitParams& initParams)
        {
        }
        void StreamLink(LbmParameters const*, geometry::FieldData&,
                        geometry::Site<geometry::FieldData> const&,
                        VarsType&, Direction) {}

        void PostStepLink(geometry::FieldData&,geometry::Site<geometry::FieldData> const&,
                          Direction) {}
    };
}

#endif
