#include "xdrReader.h"

XdrReader::XdrReader(FILE* xdrFile)
{
  xdrstdio_create (&myXdr, xdrFile, XDR_DECODE);
}

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

XdrReader::~XdrReader()
{
  xdr_destroy (&myXdr);
}