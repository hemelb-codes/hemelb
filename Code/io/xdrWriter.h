#ifndef __xdrWriter_h_
#define __xdrWriter_h_

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "writer.h"

class XdrWriter : public Writer {
 public:  
  // Method to get the current position of writing in the stream.
  unsigned int getCurrentStreamPosition() const;
  
  // Methods for formatting control
  void writeFieldSeparator();
  void writeRecordSeparator();
  
 protected:
  XDR *myXdr;
  
  // Methods to write basic types to the Xdr object.
  void _write(int const& intToWrite);
  void _write(double const& doubleToWrite);
  void _write(short const& shortToWrite);
  void _write(float const& floatToWrite);
  void _write(unsigned int const& uIntToWrite);

};

#endif //__xdrWriter_h_
