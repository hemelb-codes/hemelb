#include "xdrReader.h"

// Constructor to create an Xdr object based on a file.
XdrReader::XdrReader(FILE* xdrFile)
{
  xdrstdio_create (&myXdr, xdrFile, XDR_DECODE);
}

// Functions to read out the next bit of the file as a certain type.
void XdrReader::readDouble(double& outDouble)
{
  xdr_double(&myXdr, &outDouble);
}

void XdrReader::readInt(int& outInt)
{
  xdr_int(&myXdr, &outInt);
}

void XdrReader::readUnsignedInt(unsigned int& outUInt)
{
  xdr_u_int (&myXdr, &outUInt);
}

// Destructor to get rid of any resources used by the Xdr object. This class doesn't create
// the file object, so it doesn't free it either.
XdrReader::~XdrReader()
{
  xdr_destroy (&myXdr);
}
