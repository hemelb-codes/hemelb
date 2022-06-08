// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H

#include <cstdio>
#include <memory>
#include <vector>
#include "io/writers/xdr/XdrReader.h"

namespace hemelb::io::writers::xdr
{

  // Deserialise from a file given by path.
  class XdrFileReader : public XdrReader {
  public:
    XdrFileReader(const std::string& fn);
    ~XdrFileReader() override = default;
    unsigned GetPosition() override;
  protected:
    const char* get_bytes(size_t n) override;
  private:
    using deleter_func_t = void(*)(std::FILE*);
    std::unique_ptr<std::FILE, deleter_func_t> fh;
    // Buffer for holding read data
    std::vector<char> buf;
  };

}

#endif  // HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H
