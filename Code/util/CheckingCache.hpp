// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_CHECKINGCACHE_HPP
#define HEMELB_UTIL_CHECKINGCACHE_HPP

#include "log/Logger.h"
#include "util/CheckingCache.h"
#include "util/Cache.hpp"

namespace hemelb
{
  namespace util
  {
    template<typename CacheType>
    CheckingCache<CacheType>::CheckingCache(const lb::SimulationState& simulationState, unsigned long size) :
        Cache<CacheType>(size), simulationState(simulationState), lastUpdate(size, 0)
    {

    }

    template<typename CacheType>
    const CacheType& CheckingCache<CacheType>::Get(unsigned long index) const
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        if (lastUpdate[index] != simulationState.GetTimeStep())
        {
          log::Logger::Log<log::Warning, log::OnePerCore>("The cache was out of date at index %i", index);
        }
      }

      return Cache<CacheType>::Get(index);
    }

    template<typename CacheType>
    void CheckingCache<CacheType>::Put(unsigned long index, const CacheType& item)
    {
      lastUpdate[index] = simulationState.GetTimeStep();
      Cache<CacheType>::Put(index, item);
    }

    template<typename CacheType>
    void CheckingCache<CacheType>::Reserve(unsigned long size)
    {
      lastUpdate.resize(size, 0);
      Cache<CacheType>::Reserve(size);
    }

  }
}

#endif /* HEMELB_UTIL_CHECKINGCACHE_HPP */
