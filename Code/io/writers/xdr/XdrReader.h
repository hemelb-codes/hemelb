// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRREADER_H

#include <cstdint>
#include <string>

#include "Exception.h"
#include "io/writers/xdr/XdrSerialisation.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
        // Class to read XDR data
        class XdrReader
        {
	public:
	  // destructor.
	  virtual ~XdrReader() {
	  }

	  // Functions for reading the next bit of the stream.
	  template<class T>
	  bool read(T& val) {
	    constexpr auto n = detail::xdr_serialised_size<T>();
	    auto buf = get_bytes(n);
	    detail::xdr_deserialise(val, buf);
	    return true;
	  }

	  template <class T>
	  T read() {
	    T ans;
	    if (!read<T>(ans)) {
	      throw Exception() << "Error reading type from XDR";
	    }
	    return ans;
	  }

	  // Get the position in the stream.
	  virtual unsigned GetPosition() = 0;

	protected:
	  // Get some bytes from the underlying storage
	  virtual const char* get_bytes(size_t n) = 0;
        };

	template<>
	inline bool XdrReader::read(std::string& val) {
	  uint32_t len = 0;
	  read(len);
	  // p == padded
	  uint32_t plen = 4 * ((len - 1)/4 + 1);
	  auto pstr = get_bytes(plen);
	  val.assign(pstr, len);
	  return true;
	}

      } // namespace xdr
    } // namespace writers
  }
}

#endif // HEMELB_IO_WRITERS_XDR_XDRREADER_H
