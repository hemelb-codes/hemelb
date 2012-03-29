#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       const geometry::LatticeData& latticeData) :
      densityCache(simState, latticeData.GetLocalFluidSiteCount()),
          velocityCache(simState, latticeData.GetLocalFluidSiteCount()),
          shearStressCache(simState, latticeData.GetLocalFluidSiteCount()),
          vonMisesStressCache(simState, latticeData.GetLocalFluidSiteCount()),
          shearRateCache(simState, latticeData.GetLocalFluidSiteCount()),
          simulationState(simState),
          siteCount(latticeData.GetLocalFluidSiteCount())
    {
      ResetRequirements();
    }

    void MacroscopicPropertyCache::ResetRequirements()
    {
      densityCache.UnsetRefreshFlag();
      velocityCache.UnsetRefreshFlag();
      vonMisesStressCache.UnsetRefreshFlag();
      shearStressCache.UnsetRefreshFlag();
      shearRateCache.UnsetRefreshFlag();
    }

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
  }
}

