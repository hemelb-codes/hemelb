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
#include <iostream>
#include <tuple>

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
            iterator(ITERATOR const &iter, size_t index) : index(index), iter(iter)
            {
            }
            iterator & operator++()
            {
              ++index;
              ++iter;
              return *this;
            }
            iterator operator++(int)
            {
              iterator const result(iter.iter);
              iter++;
              return result;
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
    template<class CONTAINER> auto enumerate(CONTAINER &&x) -> Enumerate<decltype(begin(x))>
    {
      return {begin(x), end(x)};
    }

    //! Ranged-for loops with enumeration
    template<class CONTAINER> auto cenumerate(CONTAINER const &x) -> Enumerate<decltype(begin(x))>
    {
      return {begin(x), end(x)};
    }

    //! For-ranged loop over several containers in parallel
    template<class ... CONTAINER> class Zip
    {
      public:
        //! Iterates over all objects that are being zipped
        template<class ... ITERATOR> class iterator
        {
          public:
            typedef std::tuple<typename std::iterator_traits<ITERATOR>::reference...> reference;
            //! Inherits from all iterator types
            iterator(ITERATOR const & ... args) : iter(args...)
            {
            }
            //! Inherits from all iterator types
            iterator(ITERATOR && ... args) : iter(std::move(args)...)
            {
            }
            iterator& operator++()
            {
              increment<0>();
              return *this;
            }
            reference operator*() const
            {
              return deref<sizeof...(ITERATOR) - 1, sizeof...(ITERATOR) - 1>();
            }

            //! True if any of the iterator is equal to the corresponding on in b
            bool operator==(iterator const &b) const
            {
              return IsEqual<0>(b);
            }
            bool operator!=(iterator const &b) const
            {
              return not operator==(b);
            }

          protected:
            //! Increments by recurrence
            template<size_t N>
              typename std::enable_if<N < sizeof...(ITERATOR)>::type increment()
              {
                ++std::get<N>(iter);
                increment<N+1>();
              }
            //! Final for increment by recurrence
            template<size_t N>
              typename std::enable_if<N >= sizeof...(ITERATOR)>::type increment()
              {
              }
            //! Comparison by recurrence
            template<size_t N>
              typename std::enable_if<N + 1 < sizeof...(ITERATOR), bool>::type
              IsEqual(iterator const &b) const
              {
                return std::get<N>(iter) == std::get<N>(b.iter) or IsEqual<N+1>(b);
              }
            //! Final for comparison by recurrence
            template<size_t N>
              typename std::enable_if<N + 1 == sizeof...(ITERATOR), bool>::type
              IsEqual(iterator const &b) const
              {
                return std::get<N>(iter) == std::get<N>(b.iter);
              }
            //! Comparison by recurrence
            template<size_t i, size_t ... N>
              typename std::enable_if<i == 0, reference>::type deref() const
              {
                return std::tie(*std::get<N>(iter)...);
              }
            //! Comparison by recurrence
            template<size_t i, size_t ... N>
              typename std::enable_if<i != 0, reference>::type deref() const
              {
                return deref<i - 1, i - 1, N...>();
              }

            //! Holds each iterator
            std::tuple<ITERATOR...> iter;
        };

        Zip(CONTAINER && ... cont) : first(cont.begin()...), last(cont.end()...)
        {
        };
        iterator<decltype(std::declval<CONTAINER>().begin())...> begin() const
        {
          return first;
        }
        iterator<decltype(std::declval<CONTAINER>().end())...> end() const
        {
          return last;
        }
      protected:
        iterator<decltype(std::declval<CONTAINER>().begin())...> const first;
        iterator<decltype(std::declval<CONTAINER>().begin())...> const last;
    };

    //! \brief Ranged-for loops over sets of containers in parallel
    //! \code
    //! for(auto const && item: zip(container0, container2, containerN))
    //! {
    //!   std::get<0>(item);
    //!   std::get<1>(item);
    //! }
    //! \endcode
    template<class ... CONTAINER> auto zip(CONTAINER && ... x) -> Zip<CONTAINER...>
    {
      return {x...};
    }
    //! \brief Ranged-for loops over sets of (const) containers in parallel
    //! \code
    //! for(auto const && item: zip(container0, container2, containerN))
    //! {
    //!   std::get<0>(item);
    //!   std::get<1>(item);
    //! }
    //! \endcode
    template<class ... CONTAINER> auto czip(CONTAINER const & ... x) -> Zip<CONTAINER const &...>
    {
      return {x...};
    }

  }
}

#endif /* HEMELB_UTIL_CACHE_H */
