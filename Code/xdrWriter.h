#ifndef __xdrWriter_h_
#define __xdrWriter_h_

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "rt.h"
#include "xdrWriter.h"

class XdrWriter
{
  protected:
    // Constructor and member variable. This class should never be instantiated directly, so the 
    // constructor is only available to subclasses.
    XDR myXdr;
    XdrWriter();

  public:
    // Functions to write basix types to the Xdr object.
    void writePixel (ColPixel *col_pixel_p, void (*ColourPalette) (float value, float col[]));
    void writeInt(int* intToWrite);
    void writeDouble(double* doubleToWrite);
    void writeShort(short* shortToWrite);
    void writeFloat(float* floatToWrite);
    void writeUnsignedInt(unsigned int* unsignedIntToWrite);    

    // Function to get the current position of writing in the stream.
    unsigned int getCurrentStreamPosition();
};

#endif //__xdrWriter_h_
