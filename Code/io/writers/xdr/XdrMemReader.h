#ifndef HEMELB_IO_XDRFILEREADER_H
#define HEMELB_IO_XDRFILEREADER_H

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

        class XdrMemReader : public XdrReader
        {
          public:
            XdrMemReader(char* dataBuffer, unsigned int dataLength);

        };
      } // namespace xdr
    } // namespace writers
  }
}

#endif /* HEMELB_IO_XDRFILEREADER_H */
