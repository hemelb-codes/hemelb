#ifndef __io_XdrMemWriter_h_
#define __io_XdrMemWriter_h_

#include "io/XdrWriter.h"

namespace io {
  // Class for writing to an Xdr object in memory. Derives from a base
  // XdrWriter which implements all the reading functions.

  class XdrMemWriter : public XdrWriter {
  public:
    // Constructor and destructor for the in-memory Xdr writer.
    XdrMemWriter(char* dataBuffer, unsigned int dataLength);
    ~XdrMemWriter();
  };
}

#endif //__io_XdrMemWriter_h_
