// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRWRITER_H
#define HEMELB_IO_WRITERS_XDR_XDRWRITER_H

#include <cassert>

#include "io/writers/Writer.h"
#include "io/writers/xdr/XdrSerialisation.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
	// Base XDR writer - mainly to give a common base class for
	// the generic ones.
	class XdrWriter : public Writer {
	protected:
	  virtual void writeFieldSeparator();
	  virtual void writeRecordSeparator();
	};

	// Main XDR writer class template
	//
	// ByteOutputIterator - the iterator we are going to use to
	// store bytes, should dereference to char (technically we
	// want:
	// 
	// *std::declval<ByteOutputIterator>() = std::declval<char>()
	// 
	// to be valid).
	//
	// Resource - a resource we might potentially use to store data
	template<typename ByteOutputIterator,
		 typename Resource = detail::Null>
	class XdrMetaWriter : public XdrWriter
	{
	protected:
	  // Let resource in, if need be
	  friend Resource;
	  Resource res;

	  ByteOutputIterator current;

	  // Some of the following could be an empty struct depending
	  // on the traits for the iterator
	  using boi_traits = detail::boi_traits<ByteOutputIterator>;
	  typename boi_traits::counter_type bytes_written;
	  typename boi_traits::start_type start;
	  typename boi_traits::end_type end;

	  // Special constructor for difficult derived types to call.
	  //
	  // F and G are invokable to initialise the Resource and
	  // current BOI. Must have the following overloades:
	  //
	  // Resource F::operator()();
	  // ByteOperatorIterator G::operator()(Resource&);
	  //
	  // This is a bit of a hack to make Xdr{Vector, File}Writer
	  // easy...
	  template <typename F, typename G>
	  XdrMetaWriter(F f, G g) :
	    res(f()),
	    current(g(res)),
	    bytes_written(0),
	    start(current),
	    end()
	  {
	  }

	public:
	  // In the case of the Resource being Null or default
	  // constructible have two nice easy constructors.
	  //
	  // We have an iterator that supports indefinite insertion
	  XdrMetaWriter(ByteOutputIterator start_) :
	    current(start_),
	    bytes_written(0),
	    start(start_),
	    end()
	  {
	  }
	  // We have a finite amount of space between start and end
	  XdrMetaWriter(ByteOutputIterator start_, ByteOutputIterator end_) :
	    current(start_),
	    bytes_written(0),
	    start(start_),
	    end(end_)
	  {
	    static_assert(detail::is_rai<ByteOutputIterator>::value,
			  "Can't use this constructor for non-RandomAccessIterators!");
	  }

	  virtual ~XdrMetaWriter() {
	  }

	  virtual unsigned int getCurrentStreamPosition() const {
	    return boi_traits::get_position(bytes_written, start, current);
	  }
	protected:
	  // Methods to simply write (no separators) which are virtual and
	  // hence must be overriden.
	  virtual void _write(int16_t const& intToWrite) {
	    write(intToWrite);
	  }
	  virtual void _write(uint16_t const& uIntToWrite) {
	    write(uIntToWrite);
	  }
	  virtual void _write(int32_t const& intToWrite) {
	    write(intToWrite);
	  }
	  virtual void _write(uint32_t const& uIntToWrite) {
	    write(uIntToWrite);
	  }
	  virtual void _write(int64_t const& intToWrite) {
	    write(intToWrite);
	  }
	  virtual void _write(uint64_t const& uIntToWrite) {
	    write(uIntToWrite);
	  }

	  virtual void _write(double const& doubleToWrite) {
	    write(doubleToWrite);
	  }
	  virtual void _write(float const& floatToWrite) {
	    write(floatToWrite);
	  }

	  virtual void _write(const std::string& stringToWrite) {
	    // The standard defines a string of n (numbered 0 through
	    // n-1) ASCII bytes to be the number n encoded as an
	    // unsigned integer (as described above), and followed by
	    // the n bytes of the string.  Byte m of the string always
	    // precedes byte m+1 of the string, and byte 0 of the
	    // string always follows the string's length.  If n is not
	    // a multiple of four, then the n bytes are followed by
	    // enough (0 to 3) residual zero bytes, r, to make the
	    // total byte count a multiple of four.
	    const uint32_t len = stringToWrite.size();
	    write(len);
	    current = std::copy(stringToWrite.begin(), stringToWrite.end(), current);
	    // Padding
	    const auto req_nwords = (len - 1)/4 + 1;
	    for (auto rem = 4 * req_nwords - len; rem; --rem)
	      *current++ = 0;
	    bytes_written += 4*req_nwords;
	  }

	  template <typename T>
	  void write(T const& valToWrite) {
	    constexpr auto buf_size = detail::xdr_serialised_size<T>();
	    // If we have an end then we must have space to store the serialised value
	    assert(boi_traits::check_space(end, current, buf_size));

	    char buf[buf_size];
	    detail::xdr_serialise(valToWrite, buf);
	    current = std::copy(buf, buf + buf_size, current);
	    bytes_written += buf_size;
	  }
	};

      } // namespace xdr
    } // namespace writers

    template <typename ItT>
    writers::xdr::XdrMetaWriter<ItT> MakeXdrWriter(ItT start) {
      return writers::xdr::XdrMetaWriter<ItT>(start);
    }
    template <typename ItT>
    writers::xdr::XdrMetaWriter<ItT> MakeXdrWriter(ItT start, ItT end) {
      return writers::xdr::XdrMetaWriter<ItT>(start, end);
    }

  }
}

#endif //HEMELB_IO_WRITERS_XDR_XDRWRITER_H
