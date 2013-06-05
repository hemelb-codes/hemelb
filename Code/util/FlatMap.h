#ifndef HEMELB_UTIL_FLATMAP_H
#define HEMELB_UTIL_FLATMAP_H
/*
 * Using boost::container::flat_map/flat_multimap can be more efficient than
 * std::map/multimap, but we don't want to be entirely dependent on boost.
 *
 * Use template typedef metafunctions to hide the implementation detail.
 */

#ifdef HEMELB_USE_BOOST
// If boost is available.
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

#else

// Boost is not available.
#include <map>
namespace hemelb
{
  namespace util
  {
    template<class Key, class T, class Compare = std::less<Key>, class Alloc = std::allocator<
        std::pair<const Key, T> > > struct FlatMap
    {
        typedef std::map<Key, T, Compare, Alloc> Type;
    };

    template<class Key, class T, class Compare = std::less<Key>, class Alloc = std::allocator<
        std::pair<const Key, T> > >
    struct FlatMultiMap
    {
        typedef std::multimap<Key, T, Compare, Alloc> Type;
    };

  }
}
#endif

#endif // HEMELB_UTIL_FLATMAP_H
