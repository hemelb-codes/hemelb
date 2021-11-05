// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRSERIALISATION_H
#define HEMELB_IO_WRITERS_XDR_XDRSERIALISATION_H

#include <arpa/inet.h>

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
	namespace detail {
	  // Implementation details!

	  // Return the XDR representation size of a type
	  template <typename T>
	  constexpr size_t xdr_serialised_size() {
	    return 4 * ((sizeof(T)-1)/4 + 1);
	  }

	  // Overload set for XDR primitives to serialise a value into a buffer
	  // Basic signature is :
	  //  void xdr_serialise(const T& value, char* buffer);

	  // Signed and unsigned integers of 32 bits or less
	  template<typename T>
	  typename std::enable_if< std::is_integral<T>::value && sizeof(T) <= 4, void >::type
	  xdr_serialise(const T& val, char* dest_buf)
	  {
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    out[0] = htonl(val);
	  }
	    
	  // Signed and unsigned 64 bit ints
	  template<typename T>
	  typename std::enable_if< std::is_integral<T>::value && sizeof(T) == 8 >::type
	  xdr_serialise(const T& val, char* dest_buf)
	  {
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    auto msw = uint32_t(val >> 32);
	    auto lsw = uint32_t(val);
	    out[0] = htonl(msw);
	    out[1] = htonl(lsw);
	  }

	  // 32 bit floats
	  inline void xdr_serialise(const float& val, char* dest_buf)
	  {
	    auto data = reinterpret_cast<const uint32_t*>(&val);
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    out[0] = htonl(*data);
	  }
	  // 64 bit doubles
	  inline void xdr_serialise(const double& val, char* dest_buf)
	  {
	    auto data = reinterpret_cast<const uint64_t*>(&val);
	    xdr_serialise(*data, dest_buf);
	  }
	  // End of xdr_serialise overloads

	  // Overload set for XDR primitives to deserialise a value from a buffer
	  // Basic signature is :
	  //  void xdr_deserialise(T& dest, const char* buffer);
	  // Signed and unsigned integers of 32 bits or less
	  template<typename T>
	  typename std::enable_if< std::is_integral<T>::value && sizeof(T) <= 4, void >::type
	  xdr_deserialise(T& val, const char* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    val = ntohl(*in);
	  }
	    
	  // Signed and unsigned 64 bit ints
	  template<typename T>
	  typename std::enable_if< std::is_integral<T>::value && sizeof(T) == 8 >::type
	  xdr_deserialise(T& val, const char* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    uint32_t msw = ntohl(in[0]);
	    uint32_t lsw = ntohl(in[1]);
	    val = T(msw) << 32 ^ lsw;
	  }

	  // 32 bit floats
	  inline void xdr_deserialise(float& val, const char* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    auto out = reinterpret_cast<uint32_t*>(&val);
	    *out = ntohl(*in);
	  }

	  // 64 bit doubles
	  inline void xdr_deserialise(double& val, const char* src_buf)
	  {
	    auto out = reinterpret_cast<uint64_t*>(&val);
	    xdr_deserialise(*out, src_buf);
	  }
	  // End of xdr_serialise overloads
	  
	  
	  // A completely empty struct to be used as a placeholder
	  struct Null {
	    template <typename... Ts>
	    constexpr Null(Ts...) {
	    }
	    inline Null& operator+=(size_t) { return *this; }
	    constexpr bool operator!() const {
	      return true;
	    }
	    constexpr char* operator*() const {
	      return nullptr;
	    }
	  };

	  // Is the parameter a random access iterator?
	  // Answer in static member value
	  template <typename I>
	  struct is_rai {
	    static constexpr bool value = std::is_same<
	      typename std::iterator_traits<I>::iterator_category,
	      std::random_access_iterator_tag
	      >::value;
	  };

	  template <typename ItT>
	  struct boi_traits_nonaddable
	  {
	    static constexpr bool can_add = false;
	    // No need to track start
	    using start_type = Null;
	    // Or end - indeed without `operator<=` it won't work
	    using end_type = Null;
	    // Instead we count the number of bytes we write
	    using counter_type = size_t;
	    // Assume that we can write forever - e.g. to a file
	    static constexpr bool check_space(const end_type& end, const ItT& current, size_t buf_size) {
	      return true;
	    }
	    // The current position is just the number of bytes written
	    static unsigned get_position(const counter_type& bytes_written,
					 const start_type& start,
					 const ItT& current) {
	      return bytes_written;
	    }
	  };

	  // Here we specialise for the case when we do have random
	  // access iterators
	  template <typename ItT>
	  struct boi_traits_addable
	  {
	    static constexpr bool can_add = true;
	    // Start is just the same (fingers crossed iterators
	    // aren't invalidated...)
	    using start_type = ItT;
	    // End is an optional iter
	    //
	    // TODO: when we move to C++17, use std::optional
	    // using end_type = std::optional<ItT>;
	    //
	    // For now, use a super simple, minimal optional (to avoid boost)
	    class end_type {
	      ItT value;
	      bool valid = false;
	    public:
	      end_type(ItT i) : value(i), valid(true) {
	      }
	      operator bool() const {
		return valid;
	      }
	      const ItT& operator*() const {
		return value;
	      }
	    };

	    // No counter needed (Null supports operator+=(size_t))
	    using counter_type = Null;
	    // Check that we have (unlimited) space
	    static bool check_space(const end_type& end, const ItT& current, size_t buf_size) {
	      return !end || (current + buf_size <= *end);
	    }
	    // The current position is how far we got from the start
	    static unsigned get_position(const counter_type& bytes_written,
					 const start_type& start,
					 const ItT& current) {
	      return current - start;
	    }
	  };

	  // Byte Output Iterator traits - to support XdrMetaWriter
	  //
	  // We have two cases: when the iterators support random
	  // access and when they don't.
	  template <typename I>
	  using boi_traits = typename std::conditional<
	    is_rai<I>::value,
	    boi_traits_addable<I>,
	    boi_traits_nonaddable<I>
	    >::type;

	} // namespace detail
      }
    }
  }
}

#endif
