
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_REFRESHABLECACHE_H
#define HEMELB_UTIL_REFRESHABLECACHE_H

#include "util/CheckingCache.h"

namespace hemelb
{
  namespace util
  {
    /**
     * A cache that includes a flag for whether or not it requires refreshing with data.
     * CacheType is the type of object being cached.
     */
    template<typename CacheType>
    class RefreshableCache : public CheckingCache<CacheType>
    {
      public:
        /**
         * Constructor, produces a cache of the given size and initialises the refresh flag
         * to false.
         * @param simulationState
         * @param size
         */
        RefreshableCache(const lb::SimulationState& simulationState, unsigned long size);

        /**
         * Set the cache to require a refresh.
         */
        void SetRefreshFlag();

        /**
         * Set the cache to not require a refresh.
         */
        void UnsetRefreshFlag();

        /**
         * True if the cache requires refreshing. False otherwise.
         * @return
         */
        bool RequiresRefresh() const;

      private:
        /**
         * Boolean to indicate whether the cache needs refreshing.
         */
        bool requiresRefreshing;
        /**
         * The size of cache that may be required.
         */
        unsigned long cacheSize;
    };
  }
}

#endif /* HEMELB_UTIL_REFRESHABLECACHE_H */
