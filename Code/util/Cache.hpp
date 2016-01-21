
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CACHE_HPP
#define HEMELB_UTIL_CACHE_HPP

namespace hemelb
{
  namespace util
  {
    template<typename CacheType>
    Cache<CacheType>::Cache(unsigned long size) :
        items(size)
    {

    }

    template<typename CacheType>
    const CacheType& Cache<CacheType>::Get(unsigned long index) const
    {
      return items[index];
    }

    template<typename CacheType>
    void Cache<CacheType>::Put(unsigned long index, const CacheType item)
    {
      items[index] = item;
    }

    template<typename CacheType>
    void Cache<CacheType>::Reserve(unsigned long size)
    {
      items.resize(size);
    }
  }
}

#endif /* HEMELB_UTIL_CACHE_HPP */
