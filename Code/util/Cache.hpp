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
      items.reserve(size);
    }
  }
}

#endif /* HEMELB_UTIL_CACHE_HPP */
