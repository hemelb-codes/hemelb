#include <rpc/types.h>
#include <rpc/xdr.h>

#include "xdrWriter.h"

// This constructor only exists so it can be set to protected to
// prevent instances of the base class being made.  The alternative of
// declaring a pure virtual constructor or destructor doesn't
// work. Factory methods might be a potential improvement here.
XdrWriter::XdrWriter() {
}

// Functions to write simple types out to the Xdr stream.
// Sadly templating is non-trivial due to XDR's naming.
void XdrWriter::_write(int const& value) {
  xdr_int(myXdr, const_cast<int *>(&value));
}

void XdrWriter::_write(double const& doubleToWrite) {
  xdr_double (myXdr, const_cast<double *>(&doubleToWrite));
}

void XdrWriter::_write(short const& shortToWrite) {
  xdr_short (myXdr, const_cast<short *>(&shortToWrite));
}

void XdrWriter::_write(float const& floatToWrite) {
  xdr_float(myXdr, const_cast<float *>(&floatToWrite));
}

void XdrWriter::_write(unsigned int const& uIntToWrite) {
  xdr_u_int (myXdr, const_cast<unsigned int *>(&uIntToWrite));
}


// Method to get the current position in the stream.
unsigned int XdrWriter::getCurrentStreamPosition() const {
  return xdr_getpos(myXdr);
}

// No field/record separators in XDR files
void XdrWriter::writeFieldSeparator() {
}
void XdrWriter::writeRecordSeparator() {
}
