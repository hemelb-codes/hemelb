#include "xdrMemWriter.h"

// Constructor for a Xdr writer held in a memory buffer.
XdrMemWriter::XdrMemWriter(char* dataBuffer, uint dataLength)
{
  xdrmem_create(&myXdr, dataBuffer, dataLength, XDR_ENCODE);
}

// Destructor for the class.
XdrMemWriter::~XdrMemWriter()
{
  xdr_destroy(&myXdr);
}
