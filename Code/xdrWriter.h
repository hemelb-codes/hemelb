#ifndef __xdrWriter_h_
#define __xdrWriter_h_

#include <rpc/types.h>
#include <rpc/xdr.h>

//TODO TEMP HACK TO MAKE ColPixel usage below valid
#include "config.h"
#include "xdrWriter.h"

class XdrWriter
{
  protected:
    XDR myXdr;
    XdrWriter();

  public:
    void writePixel (ColPixel *col_pixel_p, void (*ColourPalette) (float value, float col[]));
    void writeInt(int* intToWrite);
    void writeDouble(double* doubleToWrite);
    void writeShort(short* shortToWrite);
    void writeFloat(float* floatToWrite);
    void writeUnsignedInt(unsigned int* unsignedIntToWrite);
};

#endif //__xdrWriter_h_