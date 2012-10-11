#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       const geometry::LatticeData& latticeData) :
      densityCache(simState, latticeData.GetLocalFluidSiteCount()),
      velocityCache(simState, latticeData.GetLocalFluidSiteCount()),
      wallShearStressMagnitudeCache(simState, latticeData.GetLocalFluidSiteCount()),
      vonMisesStressCache(simState, latticeData.GetLocalFluidSiteCount()),
      shearRateCache(simState, latticeData.GetLocalFluidSiteCount()),
      stressTensorCache(simState, latticeData.GetLocalFluidSiteCount()),
      tractionVectorCache(simState, latticeData.GetLocalFluidSiteCount()),
      tangentialProjectionTractionVectorCache(simState, latticeData.GetLocalFluidSiteCount()),
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
      wallShearStressMagnitudeCache.UnsetRefreshFlag();
      shearRateCache.UnsetRefreshFlag();
      stressTensorCache.UnsetRefreshFlag();
      tractionVectorCache.UnsetRefreshFlag();
      tangentialProjectionTractionVectorCache.UnsetRefreshFlag();
    }

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
  }
}

