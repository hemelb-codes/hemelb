// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <cassert>
#include "io/writers/xdr/XdrFileReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Constructor to create an Xdr object based on a file.
        XdrFileReader::XdrFileReader(const std::string& fn) {
	  fh = std::fopen(fn.c_str(), "r");
	  assert(fh != nullptr);
        }

	XdrFileReader::~XdrFileReader() {
	  std::fclose(fh);
	}

	unsigned XdrFileReader::GetPosition() {
	  return std::ftell(fh);
	}
	const char* XdrFileReader::get_bytes(size_t n) {
	  buf.clear();
	  buf.resize(n);
	  auto nread = std::fread(buf.data(), 1, n, fh);
	  assert(nread == n);
	  return buf.data();
	}

      }
    }
  }
}
