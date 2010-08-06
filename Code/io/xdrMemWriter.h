#ifndef __io_xdrMemWriter_h_
#define __io_xdrMemWriter_h_

#include "io/xdrWriter.h"

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

#endif //__io_xdrMemWriter_h_
