// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UTIL_ITERATOR_H
#define HEMELB_UTIL_ITERATOR_H

#include <iterator>

namespace hemelb
{
  namespace util
  {
    //! For-ranged loop plus enumeration
    template<class ITERATOR> class Enumerate
    {
        //! Holds info about enumeration
        struct EnumerateItem
        {
          size_t index;
          typename ITERATOR::reference value;
        };

      public:
        class iterator
        {
          public:
            iterator(ITERATOR const &iter) : index(0), iter(iter)
            {
            }
            auto operator++() -> decltype(++std::declval<ITERATOR>())
            {
              ++index;
              return ++iter;
            }
            auto operator++(int) -> decltype(std::declval<ITERATOR>()++)
            {
              return iter++;
            }
            EnumerateItem operator*() const
            {
              return {index, *iter};
            }
            bool operator==(iterator const &b) const
            {
              return b.iter == iter;
            }
            bool operator!=(iterator const &b) const
            {
              return b.iter != iter;
            }
          protected:
            size_t index;
            ITERATOR iter;
        };

        Enumerate(ITERATOR const &first, ITERATOR const &last) : first(first), last(last)
        {
        };
        iterator begin() const
        {
          return first;
        }
        iterator end() const
        {
          return last;
        }
      protected:
        ITERATOR const first;
        ITERATOR const last;
    };

    //! Ranged-for loops with enumeration
    template<class CONTAINER> auto enumerate(CONTAINER &&x) -> Enumerate<decltype(x.begin())>
    {
      return {x.begin(), x.end()};
    }

  }
}

#endif /* HEMELB_UTIL_CACHE_H */
