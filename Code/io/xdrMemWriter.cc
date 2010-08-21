#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrMemWriter.h"

using namespace heme::io;

// Constructor for a Xdr writer held in a memory buffer.
XdrMemWriter::XdrMemWriter(char* dataBuffer, unsigned int dataLength) {
  xdrmem_create(myXdr, dataBuffer, dataLength, XDR_ENCODE);
}

// Destructor for the class.
XdrMemWriter::~XdrMemWriter()
{
  xdr_destroy(myXdr);
}
