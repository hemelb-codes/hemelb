#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrReader.h"

namespace hemelb
{
  namespace io
  {

    XdrReader::XdrReader() :
      mXdr(NULL)
    {
    }

    // Functions to read out the next bit of the file as a certain type.
    void XdrReader::readDouble(double& outDouble)
    {
      xdr_double(mXdr, &outDouble);
    }

    void XdrReader::readFloat(float& outDouble)
    {
      xdr_float(mXdr, &outDouble);
    }

    void XdrReader::readInt(int& outInt)
    {
      xdr_int(mXdr, &outInt);
    }

    void XdrReader::readUnsignedInt(unsigned int& outUInt)
    {
      xdr_u_int(mXdr, &outUInt);
    }

    unsigned int XdrReader::GetPosition()
    {
      return xdr_getpos(mXdr);
    }

    // Destructor to get rid of any resources used by the Xdr object. This class doesn't create
    // the file object, so it doesn't free it either.
    XdrReader::~XdrReader()
    {
      xdr_destroy (mXdr);
    }

  }
}
