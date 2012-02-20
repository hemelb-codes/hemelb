#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       const geometry::LatticeData& latticeData) :
        simulationState(simState), densityCache(latticeData.GetLocalFluidSiteCount()), velocityCache(latticeData.GetLocalFluidSiteCount()), stressCache(latticeData.GetLocalFluidSiteCount()), lastCacheUpdate(CacheCount,
                                                                                                                                                                                                               std::vector<
                                                                                                                                                                                                                   unsigned long>(latticeData.GetLocalFluidSiteCount(),
                                                                                                                                                                                                                                  0))
    {

    }

    void MacroscopicPropertyCache::SetDensity(site_t siteId, distribn_t value)
    {
      densityCache[siteId] = value;
      lastCacheUpdate[DensityCache][siteId] = simulationState.GetTimeStep();
    }

    void MacroscopicPropertyCache::SetVelocity(site_t siteId, distribn_t value)
    {
      velocityCache[siteId] = value;
      lastCacheUpdate[VelocityCache][siteId] = simulationState.GetTimeStep();
    }

    void MacroscopicPropertyCache::SetStress(site_t siteId, distribn_t value)
    {
      stressCache[siteId] = value;
      lastCacheUpdate[StressCache][siteId] = simulationState.GetTimeStep();
    }

    distribn_t MacroscopicPropertyCache::GetDensity(site_t siteId) const
    {
      // Note that we use the 1-indexed timestep, so it's safe to initialise the last
      // update timestep to be 0.
      if (log::Logger::ShouldDisplay<log::Warning>())
      {
        if (lastCacheUpdate[DensityCache][siteId] != simulationState.GetTimeStep())
        {
          log::Logger::Log<log::Warning, log::OnePerCore>("The density cache was out of date for site %i", siteId);
        }
      }

      return densityCache[siteId];
    }

    const util::Vector3D<distribn_t>& MacroscopicPropertyCache::GetVelocity(site_t siteId) const
    {
      // Note that we use the 1-indexed timestep, so it's safe to initialise the last
      // update timestep to be 0.
      if (log::Logger::ShouldDisplay<log::Warning>())
      {
        if (lastCacheUpdate[VelocityCache][siteId] != simulationState.GetTimeStep())
        {
          log::Logger::Log<log::Warning, log::OnePerCore>("The velocity cache was out of date for site %i", siteId);
        }
      }

      return velocityCache[siteId];
    }

    distribn_t MacroscopicPropertyCache::GetStress(site_t siteId) const
    {
      // Note that we use the 1-indexed timestep, so it's safe to initialise the last
      // update timestep to be 0.
      if (log::Logger::ShouldDisplay<log::Warning>())
      {
        if (lastCacheUpdate[StressCache][siteId] != simulationState.GetTimeStep())
        {
          log::Logger::Log<log::Warning, log::OnePerCore>("The stress cache was out of date for site %i", siteId);
        }
      }

      return stressCache[siteId];
    }
  }
}

