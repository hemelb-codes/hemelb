#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrReader.h"

namespace hemelb
{
  namespace io
  {

    XdrReader::XdrReader()
    {
    }

    // Functions to read out the next bit of the file as a certain type.
    bool XdrReader::readDouble(double& outDouble)
    {
      return xdr_double(&mXdr, &outDouble);
    }

    bool XdrReader::readFloat(float& outDouble)
    {
      return xdr_float(&mXdr, &outDouble);
    }

    bool XdrReader::readInt(int& outInt)
    {
      return xdr_int(&mXdr, &outInt);
    }

    bool XdrReader::readUnsignedInt(unsigned int& outUInt)
    {
      return xdr_u_int(&mXdr, &outUInt);
    }

    unsigned int XdrReader::GetPosition()
    {
      return xdr_getpos(&mXdr);
    }

    // Returns false on failure
    bool XdrReader::SetPosition(unsigned int iPosition)
    {
      return xdr_setpos(&mXdr, iPosition);
    }

    // Destructor to get rid of any resources used by the Xdr object. This class doesn't create
    // the file object, so it doesn't free it either.
    XdrReader::~XdrReader()
    {
      xdr_destroy (&mXdr);
    }

  }
}
