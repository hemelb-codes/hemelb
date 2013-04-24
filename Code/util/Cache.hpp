// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
