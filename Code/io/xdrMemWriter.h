#ifndef HEME_IO_XDRMEMWRITER_H
#define HEME_IO_XDRMEMWRITER_H

#include "io/XdrWriter.h"

namespace heme
{
  namespace io
  {
    // Class for writing to an Xdr object in memory. Derives from a
    // base XdrWriter which implements all the writing functions.
    
    class XdrMemWriter : public XdrWriter {
    public:
      // Constructor and destructor for the in-memory Xdr writer.
      XdrMemWriter(char* dataBuffer, unsigned int dataLength);
      ~XdrMemWriter();
    };
  }
}
#endif // HEME_IO_XDRMEMWRITER_H
