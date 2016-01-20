
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CACHE_H
#define HEMELB_UTIL_CACHE_H

#include <vector>

namespace hemelb
{
  namespace util
  {
    /**
     * Generic class for implementing a cache.
     * CacheType is the type of the cached object.
     */
    template<typename CacheType>
    class Cache
    {
      public:
        /**
         * Constructor, takes the cache size.
         * @param size
         */
        Cache(unsigned long size);
        /**
         * Gets the cached object from the case for the given index.
         * @param index
         * @return
         */
        const CacheType& Get(unsigned long index) const;
        /**
         * Inserts the given object into the cache at the given index.
         * @param index
         * @param item
         */
        void Put(unsigned long index, const CacheType item);

      protected:
        /**
         * Resizes the cache to the given size.
         * This lets us minimise the memory used by each cache until they're used.
         * @param size
         */
        void Reserve(unsigned long size);
      private:
        /**
         * The actual cached items.
         */
        std::vector<CacheType> items;
    };
  }
}

#endif /* HEMELB_UTIL_CACHE_H */
