#ifndef HEMELB_IO_XDRFILEREADER_H
#define HEMELB_IO_XDRFILEREADER_H

#include <cstdio>

#include "io/XdrReader.h"

namespace hemelb
{
  namespace io
  {

    class XdrFileReader : public XdrReader
    {
      public:
        XdrFileReader(FILE* xdrFile);
    };

  }

}

#endif /* HEMELB_IO_XDRFILEREADER_H */
