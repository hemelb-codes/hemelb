#ifndef HEMELB_IO_XDRMEMWRITER_H
#define HEMELB_IO_XDRMEMWRITER_H

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Class for writing to an Xdr object in memory. Derives from a
        // base XdrWriter which implements all the writing functions.

        class XdrMemWriter : public XdrWriter
        {
          public:
            // Constructor and destructor for the in-memory Xdr writer.
            XdrMemWriter(char* dataBuffer, unsigned int dataLength);
            ~XdrMemWriter();
        };

      } // namespace xdr
    } // namespace writers
  }
}
#endif // HEMELB_IO_XDRMEMWRITER_H
