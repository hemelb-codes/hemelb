#ifndef HEMELB_LB_MACROSCOPICPROPERTYCACHE_H
#define HEMELB_LB_MACROSCOPICPROPERTYCACHE_H

#include <vector>
#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "units.h"

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
         * Set the cached value for density at a given site.
         * @param siteId A contiguous site id on this core.
         * @param value The density there.
         */
        void SetDensity(site_t siteId, distribn_t value);

        /**
         * Set the cached value for velocity at a given site.
         * @param siteId A contiguous site id on this core.
         * @param value The velocity there.
         */
        void SetVelocity(site_t siteId, distribn_t value);

        /**
         * Set the cached value for stress at a given site.
         * @param siteId A contiguous site id on this core.
         * @param value The stress there.
         */
        void SetStress(site_t siteId, distribn_t value);

        /**
         * Returns the cached value of the density at a local contiguous site id
         * (and will check whether the cache is out of date at certain logging levels).
         * @param siteId The local contiguous site id.
         * @return The density there.
         */
        distribn_t GetDensity(site_t siteId) const;

        /**
         * Returns the cached value of the velocity at a local contiguous site id
         * (and will check whether the cache is out of date at certain logging levels).
         * @param siteId The local contiguous site id.
         * @return The velocity there.
         */
        const util::Vector3D<distribn_t>& GetVelocity(site_t siteId) const;

        /**
         * Returns the cached value of the stress at a local contiguous site id
         * (and will check whether the cache is out of date at certain logging levels).
         * @param siteId The local contiguous site id.
         * @return The stress there.
         */
        distribn_t GetStress(site_t siteId) const;

      private:
        /**
         * Enumerates each of the different caches we have. This is primarily for accessing
         * the lastCacheUpdate value.
         */
        enum
        {
          DensityCache = 0, //!< DensityCache
          VelocityCache = 1,//!< VelocityCache
          StressCache = 2, //!< StressCache

          CacheCount
        //!< CacheCount
        } CacheType;

        /**
         * The state of the simulation, including the number of timesteps passed.
         */
        const SimulationState& simulationState;
        /**
         * The cache of densities for each fluid site on this core.
         */
        std::vector<distribn_t> densityCache;
        /**
         * The cache of velocities for each fluid site on this core.
         */
        std::vector<util::Vector3D<distribn_t> > velocityCache;
        /**
         * The cache of stresses for each fluid site on this core.
         */
        std::vector<distribn_t> stressCache;
        /**
         * The timesteps at which each of the caches was updated on each local fluid site.
         */
        std::vector<std::vector<unsigned long> > lastCacheUpdate;
    };
  }
}

#endif /* HEMELB_LB_MACROSCOPICPROPERTYCACHE_H */
