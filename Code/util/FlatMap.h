
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_UTIL_FLATMAP_H
#define HEMELB_UTIL_FLATMAP_H

#include <boost/container/flat_map.hpp>

namespace hemelb
{
  namespace util
  {
    template<class Key, class T, class Compare = std::less<Key>, class Alloc = std::allocator<
    std::pair<Key, T> > >
    struct FlatMap
    {
      typedef boost::container::flat_map<Key, T, Compare, Alloc> Type;
    };

    template<class Key, class T, class Compare = std::less<Key>, class Alloc = std::allocator<
    std::pair<Key, T> > >
    struct FlatMultiMap
    {
      typedef boost::container::flat_multimap<Key, T, Compare, Alloc> Type;
    };
  }
}

#endif // HEMELB_UTIL_FLATMAP_H
