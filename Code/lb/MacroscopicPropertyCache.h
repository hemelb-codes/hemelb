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
        MacroscopicPropertyCache(const SimulationState& simState, const geometry::LatticeData& latticeData);

        void SetDensity(site_t siteId, distribn_t value);
        void SetVelocity(site_t siteId, distribn_t value);
        void SetStress(site_t siteId, distribn_t value);

        distribn_t GetDensity(site_t siteId) const;
        const util::Vector3D<distribn_t>& GetVelocity(site_t siteId) const;
        distribn_t GetStress(site_t siteId) const;

      private:
        enum
        {
          DensityCache = 0,
          VelocityCache = 1,
          StressCache = 2,

          CacheCount
        } CacheType;

        const SimulationState& simulationState;
        std::vector<distribn_t> densityCache;
        std::vector<util::Vector3D<distribn_t> > velocityCache;
        std::vector<distribn_t> stressCache;
        std::vector<std::vector<unsigned long> > lastCacheUpdate;
    };
  }
}

#endif /* HEMELB_LB_MACROSCOPICPROPERTYCACHE_H */
