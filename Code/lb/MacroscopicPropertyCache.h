// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LB_MACROSCOPICPROPERTYCACHE_H
#define HEMELB_LB_MACROSCOPICPROPERTYCACHE_H

#include <vector>
#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "units.h"
#include "util/RefreshableCache.hpp"

namespace hemelb
{
  namespace lb
  {
    class MacroscopicPropertyCache
    {
      public:
        /**
         * Constructor, the only way to create this object.
         * @param simState The simulation state, so that the cache knows when it is out of date.
         * @param latticeData Data about the lattice, so the cache knows how big it needs to be.
         * @return
         */
        MacroscopicPropertyCache(const SimulationState& simState, const geometry::LatticeData& latticeData);

        /**
         * Reset the list of cache types required to be none of them.
         */
        void ResetRequirements();

        /**
         * Returns the number of sites cached.
         * @return
         */
        site_t GetSiteCount() const;

        /**
         * The cache of densities for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> densityCache;

        /**
         * The cache of velocities for each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<distribn_t> > velocityCache;

        /**
         * The cache of wall shear stress magnitudes for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> wallShearStressMagnitudeCache;

        /**
         * The cache of Von Mises stresses for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> vonMisesStressCache;

        /**
         * The cache of shear rates for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> shearRateCache;

        /**
         * The cache of stress tensors for each fluid site on this core.
         */
        util::RefreshableCache<util::Matrix3D> stressTensorCache;

        /**
         * The cache of traction vectors for each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<LatticeStress> > tractionVectorCache;

        /**
         * The cache of projections of the traction vectors on the tangential plane of each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<LatticeStress> > tangentialProjectionTractionVectorCache;

      private:
        /**
         * The state of the simulation, including the number of timesteps passed.
         */
        const SimulationState& simulationState;

        /**
         * The number of sites.
         */
        site_t siteCount;
    };
  }
}

#endif /* HEMELB_LB_MACROSCOPICPROPERTYCACHE_H */
