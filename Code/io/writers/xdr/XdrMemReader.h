// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H

#include <vector>

#include "io/writers/xdr/XdrReader.h"

namespace hemelb::io::writers::xdr
{

  // Deserialise from some memory.
  //
  // NB: this class does not OWN the memory you pass to its
  // constructor. Caller is responsible for ensuring its lifetime
  // exceeds that of the object.
  class XdrMemReader : public XdrReader
  {
  public:
    XdrMemReader(const char* dataBuffer, unsigned int dataLength);
    XdrMemReader(const std::vector<char>& dataVec);
    ~XdrMemReader() override = default;
    unsigned GetPosition() override;
  protected:
    const char* get_bytes(size_t n) override;
  private:
    const char* start;
    const char* current;
    size_t len;
  };
}

#endif  // HEMELB_IO_XDRFILEREADER_H
