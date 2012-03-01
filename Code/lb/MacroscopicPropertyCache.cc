#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       const geometry::LatticeData& latticeData) :
      simulationState(simState), densityCache(latticeData.GetLocalFluidSiteCount()),
          velocityCache(latticeData.GetLocalFluidSiteCount()), stressCache(latticeData.GetLocalFluidSiteCount()),
          lastCacheUpdate(CacheCount, std::vector<unsigned long>(latticeData.GetLocalFluidSiteCount(), 0)),
          siteCount(latticeData.GetLocalFluidSiteCount())
    {

    }

    void MacroscopicPropertyCache::SetDensity(site_t siteId, distribn_t value)
    {
      densityCache[siteId] = value;
      lastCacheUpdate[DensityCache][siteId] = simulationState.GetTimeStep();
    }

    void MacroscopicPropertyCache::SetVelocity(site_t siteId, const util::Vector3D<distribn_t>& value)
    {
      velocityCache[siteId] = value;
      lastCacheUpdate[VelocityCache][siteId] = simulationState.GetTimeStep();
    }

    void MacroscopicPropertyCache::SetVelocity(site_t siteId,
                                               const distribn_t& v_x,
                                               const distribn_t& v_y,
                                               const distribn_t& v_z)
    {
      velocityCache[siteId].x = v_x;
      velocityCache[siteId].y = v_y;
      velocityCache[siteId].z = v_z;
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

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
  }
}

