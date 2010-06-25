#include "xdrWriter.h"

// Class for writing to an Xdr object in memory. Derives from a base XdrWriter 
// which implements all the reading functions.
class XdrMemWriter : public XdrWriter
{
  public:
    // Constructor and destructor for the in-memory Xdr writer.
    XdrMemWriter(char* dataBuffer, uint dataLength);
    ~XdrMemWriter();
};
