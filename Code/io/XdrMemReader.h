#ifndef HEMELB_IO_XDRFILEREADER_H
#define HEMELB_IO_XDRFILEREADER_H

#include <cstdio>

#include "io/XdrReader.h"

namespace hemelb
{
  namespace io
  {

    class XdrMemReader : public XdrReader
    {
      public:
        XdrMemReader(char* dataBuffer, unsigned int dataLength);

    };

  }

}

#endif /* HEMELB_IO_XDRFILEREADER_H */
