// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include "io/writers/xdr/XdrMemReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
        // Constructor to create an Xdr object based on a memory buffer
        XdrMemReader::XdrMemReader(const char* buf, unsigned int dataLength)
	  : start(buf), current(buf), len(dataLength)
        {
        }

	XdrMemReader::XdrMemReader(const std::vector<char>& dataVec)
	  : start(dataVec.data()), current(start), len(dataVec.size())
	{
	}

	XdrMemReader::~XdrMemReader() {
	}

	unsigned XdrMemReader::GetPosition() {
	  return current - start;
	}

	const char* XdrMemReader::get_bytes(size_t n) {
	  assert(GetPosition() + n <= len);
	  auto ans = current;
	  current += n;
	  return ans;
	}

      } // namespace xdr
    } // namespace writers
  }
}
