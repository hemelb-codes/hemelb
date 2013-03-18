// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_REFRESHABLECACHE_HPP
#define HEMELB_UTIL_REFRESHABLECACHE_HPP

#include "util/RefreshableCache.h"
#include "util/CheckingCache.hpp"

namespace hemelb
{
  namespace util
  {
    /**
     * NOTE: We initialise the checking cache to size 0, then expand it if needed.
     * @param simulationState
     * @param size
     */
    template<typename CacheType>
    RefreshableCache<CacheType>::RefreshableCache(const lb::SimulationState& simulationState, unsigned long size) :
        CheckingCache<CacheType>(simulationState, size), requiresRefreshing(false), cacheSize(size)
    {

    }

    template<typename CacheType>
    void RefreshableCache<CacheType>::SetRefreshFlag()
    {
      CheckingCache<CacheType>::Reserve(cacheSize);
      requiresRefreshing = true;
    }

    template<typename CacheType>
    void RefreshableCache<CacheType>::UnsetRefreshFlag()
    {
      requiresRefreshing = false;
    }

    template<typename CacheType>
    bool RefreshableCache<CacheType>::RequiresRefresh() const
    {
      return requiresRefreshing;
    }
  }
}

#endif /* HEMELB_UTIL_REFRESHABLECACHE_HPP */
