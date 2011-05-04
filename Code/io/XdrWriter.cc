#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    // Functions to write simple types out to the Xdr stream.
    // Sadly templating is non-trivial due to XDR's naming.
    void XdrWriter::_write(int const& value)
    {
      xdr_int(&mXdr, const_cast<int *> (&value));
    }

    void XdrWriter::_write(double const& doubleToWrite)
    {
      xdr_double(&mXdr, const_cast<double *> (&doubleToWrite));
    }

    void XdrWriter::_write(short const& shortToWrite)
    {
      xdr_short(&mXdr, const_cast<short *> (&shortToWrite));
    }

    void XdrWriter::_write(float const& floatToWrite)
    {
      xdr_float(&mXdr, const_cast<float *> (&floatToWrite));
    }

    void XdrWriter::_write(unsigned int const& uIntToWrite)
    {
      xdr_u_int(&mXdr, const_cast<unsigned int *> (&uIntToWrite));
    }

    void XdrWriter::_write(long const& longToWrite)
    {
      xdr_long(&mXdr, const_cast<long*> (&longToWrite));
    }

    // Method to get the current position in the stream.
    unsigned int XdrWriter::getCurrentStreamPosition() const
    {
      // The XDR function does not modify the instance, but is not const in the declaration.
      return xdr_getpos(const_cast<XDR*>(&mXdr));
    }

    // No field/record separators in XDR files
    void XdrWriter::writeFieldSeparator()
    {
    }
    void XdrWriter::writeRecordSeparator()
    {
    }

  }
}
