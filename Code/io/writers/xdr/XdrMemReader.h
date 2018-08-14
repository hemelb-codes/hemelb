
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H

#include <vector>

#include "io/writers/xdr/XdrReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        class XdrMemReader : public XdrReader
        {
          public:
            XdrMemReader(const char* dataBuffer, unsigned int dataLength);
            XdrMemReader(const std::vector<char>& dataVec);

        };
      } // namespace xdr
    } // namespace writers
  }
}

#endif /* HEMELB_IO_XDRFILEREADER_H */
