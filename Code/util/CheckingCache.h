
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CHECKINGCACHE_H
#define HEMELB_UTIL_CHECKINGCACHE_H

#include "lb/SimulationState.h"
#include "util/Cache.h"

namespace hemelb
{
  namespace util
  {
    /**
     * Cache that checks whether cached values are stale.
     * CacheType is the type of the object cached.
     */
    template<typename CacheType>
    class CheckingCache : public Cache<CacheType>
    {
      public:
        /**
         * Constructor, requires the simulation state (to get timestep id) and the size of the
         * cache.
         * @param simulationState
         * @param size
         */
        CheckingCache(const lb::SimulationState& simulationState, unsigned long size);

        /**
         * Obtain an object from the cache, checking whether it's stale.
         * NOTE: It is deliberate that this does *not* use polymorphism (efficiency is to
         * important here). We instead cover the similar method in the base class.
         * @param index
         * @return
         */
        const CacheType& Get(unsigned long index) const;

        /**
         * Inserts the given object into the cache at the given index.
         * NOTE: It is deliberate that this does *not* use polymorphism (efficiency is to
         * important here). We instead cover the similar method in the base class.
         * @param index
         * @param item
         */
        void Put(unsigned long index, const CacheType& item);

      protected:
        /**
         * Reserves enough space for the cache.
         *
         * @param size
         */
        void Reserve(unsigned long size);

        /**
         * The state of the simulation, for accessing the timestep id.
         */
        const lb::SimulationState& simulationState;

        /**
         * The last update at each cache position.
         */
        std::vector<unsigned long> lastUpdate;
    };
  }
}

#endif /* HEMELB_UTIL_CHECKINGCACHE_H */
