#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrMemWriter.h"

namespace hemelb
{
  namespace io
  {

    // Constructor for a Xdr writer held in a memory buffer.
    XdrMemWriter::XdrMemWriter(char* dataBuffer, unsigned int dataLength) {
      myXdr = new XDR;
      xdrmem_create(myXdr, dataBuffer, dataLength, XDR_ENCODE);
    }

    // Destructor for the class.
    XdrMemWriter::~XdrMemWriter()
    {
      xdr_destroy(myXdr);
      delete myXdr;
    }
    
  }
}
