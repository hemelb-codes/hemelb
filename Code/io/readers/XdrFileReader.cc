// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <cassert>
#include <cstring>

#include "io/readers/XdrFileReader.h"

namespace hemelb::io
{
  namespace {
    void fclose_deleter(std::FILE* fh) {
      std::fclose(fh);
    }

  }

  // Constructor to create an Xdr object based on a file.
  XdrFileReader::XdrFileReader(const std::string& fn) :
    fh{std::fopen(fn.c_str(), "r"), fclose_deleter} {
    if (!fh) {
      throw Exception() << "Error opening file '" << fn << "' with reason: " << std::strerror(errno);
    }
  }

  unsigned XdrFileReader::GetPosition() {
    return std::ftell(fh.get());
  }

  const char* XdrFileReader::get_bytes(size_t n) {
    buf.clear();
    buf.resize(n);
    auto nread = std::fread(buf.data(), 1, n, fh.get());
    assert(nread == n);
    return buf.data();
  }

}
