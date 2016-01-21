
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H

#include <cstdio>

#include "io/writers/xdr/XdrReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {

      namespace xdr
      {

        class XdrFileReader : public XdrReader
        {
          public:
            XdrFileReader(FILE* xdrFile);
        };

      } // namespace xdr

    } // namespace writers

  }

}

#endif /* HEMELB_IO_WRITERS_XDR_XDRFILEREADER_H */
