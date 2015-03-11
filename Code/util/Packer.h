//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_UTIL_PACKER
#define HEMELB_UTIL_PACKER

#include <vector>
#include <utility>
#include <map>
#include <type_traits>
#include <cassert>
#include <string>
#include "util/Vector3D.h"

namespace hemelb
{
  namespace util
  {
    // Advanced declaration for friend functions
    class Packer;

    //! packs object to packer
    template<class T> typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
      operator<<(Packer& buffer, T const & packme);
    //! unpacks object from packer
    template<class T> typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
      operator>>(Packer& buffer, T & packme);

    //! \brief Packs and unpacks to an internal stream
    //! \details Stores information into a buffer. The operators ``<<`` and ``>>`` to and from the
    //! packer object allow storing information. This is a fairly simple implementation. There is no
    //! support for pointers and such (eg recreating two pointers to the same object).
    class Packer
    {
      // These two should be sufficient to build other streaming operators.
      template<class T>
        friend typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
        operator<<(Packer& buffer, T const & packme);
      template<class T>
        friend typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
        operator>>(Packer& buffer, T & packme);

      public:
        //! type underlying the packer
        typedef int8_t StreamByte;
        //! type underlying the packer
        typedef std::vector<StreamByte> Buffer;
        //! type for sizes
        typedef Buffer::size_type size_type;

        //! Convenience class to sum up sizes
        class Sizer;

        //! New empty buffer
        Packer() : buffer(new Buffer), pos(buffer->begin())
        {
        }
        Packer(Buffer const &b) : buffer(new Buffer(b)), pos(buffer->begin())
        {
        }
        Packer(Buffer &&b) : buffer(new Buffer(std::move(b))), pos(buffer->begin())
        {
        }
        //! Shallow copy constructor
        Packer(Packer const &c) : buffer(c.buffer), pos(c.pos)
        {
        }

        //! pointer to data
        StreamByte* data()
        {
          return buffer->data();
        }
        //! pointer to data
        StreamByte const* data() const
        {
          return buffer->data();
        }
        //! current buffer size
        Buffer::size_type size()
        {
          return buffer->size();
        }
        //! current buffer size in bytes
        Buffer::size_type messageSize()
        {
          return buffer->size() * sizeof(StreamByte);
        }
        //! Clears current buffer
        void clear()
        {
          buffer->clear();
        }

        //! current read position at beginning
        void reset_read()
        {
          pos = buffer->begin();
        }

      protected:
        //! Object holding the packer
        std::shared_ptr<Buffer> buffer;
        //! Object holding the packer
        Buffer::const_iterator pos;
    };

    //! packs object to packer
    template<class T> typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
      operator<<(Packer& packer, T const & packme)
      {
        if(sizeof(T) == sizeof(Packer::StreamByte))
        {
          packer.buffer->push_back(*reinterpret_cast<Packer::StreamByte const*>(&packme));
        }
        else if(sizeof(T) < sizeof(Packer::StreamByte))
        {
          packer.buffer->push_back(0);
          *reinterpret_cast<T*>(&packer.buffer->back()) = packme;
        }
        else
        {
          assert(sizeof(T) % sizeof(Packer::StreamByte) == 0);
          auto const first = reinterpret_cast<Packer::StreamByte const*>(&packme);
          auto const end = first + sizeof(T) / sizeof(Packer::StreamByte);
          packer.buffer->insert(packer.buffer->end(), first, end);
        }
        return packer;
      }
    //! unpacks object from packer
    template<class T> typename std::enable_if<std::is_arithmetic<T>::value, Packer&> :: type
      operator>>(Packer& packer, T & packme)
      {
        assert(packer.pos != packer.buffer->end());
        packme = *reinterpret_cast<T const *>(&(*packer.pos));
        if(sizeof(T) <= sizeof(Packer::StreamByte))
        {
          ++packer.pos;
        }
        else
        {
          assert(sizeof(T) % sizeof(Packer::StreamByte) == 0);
          assert(packer.pos + sizeof(T) / sizeof(Packer::StreamByte) <= packer.buffer->end());
          packer.pos += sizeof(T) / sizeof(Packer::StreamByte);
        }
        return packer;
      }


    //! packing for 3d vectors
    template<class T> Packer& operator<<(Packer &packer, Vector3D<T> const &vector)
    {
      return packer << vector.x << vector.y << vector.z;
    }
    //! unpacking for 3d vectors
    template<class T> Packer& operator>>(Packer &packer, Vector3D<T> & vector)
    {
      return packer >> vector.x >> vector.y >> vector.z;
    }

    template<class T, class D>
      typename std::enable_if<std::is_default_constructible<T>::value, Packer&>::type
      operator<<(Packer &packer, std::vector<T, D> const &vector)
    {
      packer << vector.size();
      for(auto const & element: vector)
      {
        packer << element;
      }
      return packer;
    }
    template<class T, class D>
      typename std::enable_if<std::is_default_constructible<T>::value, Packer&>::type
      operator>>(Packer &packer, std::vector<T, D> &vector)
    {
      typename std::vector<T, D>::size_type n;
      packer >> n;
      vector.resize(n);
      for(auto & element: vector)
      {
        packer >> element;
      }
      return packer;
    }

    template<class A, class B>
      typename std::enable_if<
        std::is_default_constructible<A>::value
        and std::is_default_constructible<B>::value,
        Packer&
      >::type operator<<(Packer &packer, std::pair<A, B> const &pair)
    {
      return packer << pair.first << pair.second;
    }
    template<class A, class B>
      typename std::enable_if<
        std::is_default_constructible<A>::value
        and std::is_default_constructible<B>::value,
        Packer&
      >::type operator>>(Packer &packer, std::pair<A, B> &pair)
    {
      return packer >> pair.first >> pair.second;
    }

    template<class KEY, class T, class COMPARE, class ALLOC>
      typename std::enable_if<
        std::is_default_constructible<KEY>::value
        and std::is_default_constructible<T>::value,
        Packer&
      >::type operator<<(Packer &packer, std::map<KEY, T, COMPARE, ALLOC> const &map)
    {
      packer << map.size();
      for(auto const & element: map)
      {
        packer << element;
      }
      return packer;
    }
    template<class KEY, class T, class COMPARE, class ALLOC>
      typename std::enable_if<
        std::is_default_constructible<KEY>::value
        and std::is_default_constructible<T>::value,
        Packer&
      >::type operator>>(Packer &packer, std::map<KEY, T, COMPARE, ALLOC> &map)
    {
      typedef std::map<KEY, T, COMPARE, ALLOC> Map;
      typedef typename std::remove_const<typename Map::key_type>::type key_type;
      typedef typename std::remove_const<typename Map::mapped_type>::type mapped_type;
      typename Map::size_type n;
      packer >> n;
      map.clear();
      for(typename Map::size_type i(0); i < n; ++i)
      {
        std::pair<key_type, mapped_type> value;
        packer >> value;
        map.emplace(std::move(value));
      }
      return packer;
    }

    template<class T>
      Packer& operator<<(Packer& packer, std::basic_string<T> const &string)
      {
        packer << string.size();
        for(auto const & character: string)
        {
          packer << character;
        }
        return packer;
      }
    template<class T>
      Packer& operator>>(Packer& packer, std::basic_string<T> &string)
      {
        decltype(string.size()) n;
        packer >> n;
        for(decltype(n) i(0); i < n; ++i)
        {
          T character;
          packer >> character;
          string += character;
        }
        return packer;
      }

    class Packer::Sizer
    {
      public:
        Sizer(Packer::Buffer::size_type i = 0) : accumulated(i)
        {
        }
        Packer::Buffer::size_type cast() const
        {
          return accumulated;
        }
        void operator+=(Packer::Buffer::size_type i)
        {
          accumulated += i;
        }
        // Adds one size n times
        template<class T> void addNElements(Packer::Buffer::size_type n, T const& value)
        {
          if(n != 0)
          {
            Sizer other;
            other << value;
            accumulated += n * other.cast();
          }
        }
      private:
        Packer::Buffer::size_type accumulated;
    };

    namespace details
    {
      //! True if a type has fixed memory footprint
      template<class T, class V = void> struct FixedMemoryFootPrint : std::false_type
      {
      };

      //! True if a type has fixed memory footprint
      template<class T>
        struct FixedMemoryFootPrint<
          T,
          typename std::enable_if<std::is_arithmetic<T>::value, void>::type
        > : std::true_type
      {
      };

      //! True if a type has fixed memory footprint
      template<class T>
        struct FixedMemoryFootPrint<
          util::Vector3D<T>,
          typename std::enable_if<FixedMemoryFootPrint<T>::value, void>::type
        > : std::true_type
        {
        };
    }

    //! figures out size of message in bytes
    template<class T>
      typename std::enable_if<std::is_arithmetic<T>::value, Packer::Sizer&> :: type
      operator<<(Packer::Sizer& sizer, T const & packme)
      {
        sizer += std::max(sizeof(T), sizeof(Packer::StreamByte));
        return sizer;
      }
    template<class T> Packer::Sizer& operator<<(Packer::Sizer& sizer, Vector3D<T> const &v)
    {
      sizer.addNElements(3, v.x);
      return sizer;
    }
    template<class T, class D>
      typename std::enable_if<
        std::is_default_constructible<T>::value
        and details::FixedMemoryFootPrint<T>::value, 
        Packer::Sizer&
      >::type operator<<(Packer::Sizer& sizer, std::vector<T, D> const &vector)
    {
      auto const N = vector.size();
      sizer << N;
      if(N != 0)
        sizer.addNElements(N, vector.front());
      return sizer;
    }
    template<class T, class D>
      typename std::enable_if<
        std::is_default_constructible<T>::value
        and not details::FixedMemoryFootPrint<T>::value,
        Packer::Sizer&
      >::type
      operator<<(Packer::Sizer& sizer, std::vector<T, D> const &vector)
    {
      auto const N = vector.size();
      sizer << N;
      for(auto const& element: vector)
      {
        sizer << element;
      }
      return sizer;
    }
    template<class A, class B>
      typename std::enable_if<
        std::is_default_constructible<A>::value
        and std::is_default_constructible<B>::value,
        Packer::Sizer&
      >::type operator<<(Packer::Sizer& sizer, std::pair<A, B> const &pair)
    {
      return sizer << pair.first << pair.second;
    }
    template<class KEY, class T, class COMPARE, class ALLOC>
      typename std::enable_if<
        std::is_default_constructible<KEY>::value
        and std::is_default_constructible<T>::value
        and details::FixedMemoryFootPrint<T>::value
        and details::FixedMemoryFootPrint<KEY>::value,
        Packer::Sizer&
      >::type operator<<(Packer::Sizer& sizer, std::map<KEY, T, COMPARE, ALLOC> const &map)
    {
      typedef typename std::map<KEY, T, COMPARE, ALLOC>::value_type value_type;
      auto const N = map.size();
      sizer << N;
      if(N != 0)
        sizer.addNElements(N, *map.begin());
      return sizer;
    }
    template<class KEY, class T, class COMPARE, class ALLOC>
      typename std::enable_if<
        std::is_default_constructible<KEY>::value
        and std::is_default_constructible<T>::value
        and not (
            details::FixedMemoryFootPrint<T>::value and details::FixedMemoryFootPrint<KEY>::value),
        Packer::Sizer&
      >::type operator<<(Packer::Sizer& sizer, std::map<KEY, T, COMPARE, ALLOC> const &map)
    {
      typedef typename std::map<KEY, T, COMPARE, ALLOC>::value_type value_type;
      auto const N = map.size();
      sizer << N;
      for(auto const & element: map)
      {
        sizer << element;
      }
      return sizer;
    }
    template<class T>
      Packer::Sizer& operator<<(Packer::Sizer& sizer, std::basic_string<T> const &string)
      {
        auto const N = string.size();
        sizer << N;
        if(N != 0)
          sizer.addNElements(N, string[0]);
        return sizer;
      }

    namespace details
    {
      template<class T, class... Arguments>
        Packer::size_type
        packerMessageSizeImpl(Packer::Sizer &sizer, T const &value, Arguments const &... args)
      {
        sizer << value;
        packerMessageSizeImpl(sizer, args...);
        return sizer.cast();
      }
      template<class T>
        Packer::size_type packerMessageSizeImpl(Packer::Sizer &sizer, T const &value)
      {
        sizer << value;
        return sizer.cast();
      }
    }

    // Convenience function to compute message size of multiple argument
    template<class T, class... Arguments>
      typename std::enable_if<
        not std::is_same<typename std::remove_all_extents<T>::type, Packer::Sizer>::value,
        Packer::size_type
      >::type packerMessageSize(T const &value, Arguments const &... args)
    {
      Packer::Sizer sizer;
      sizer << value;
      details::packerMessageSizeImpl(sizer, args...);
      return sizer.cast();
    }

    // Convenience function to compute message size of single argument
    template<class T>
      typename std::enable_if<
        not std::is_same<typename std::remove_all_extents<T>::type, Packer::Sizer>::value,
        Packer::size_type
      >::type packerMessageSize(T const &value)
    {
      Packer::Sizer sizer;
      sizer << value;
      return sizer.cast();
    }

  } // util
} // hemelb

#endif
