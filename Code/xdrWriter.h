#ifndef __xdrWriter_h_
#define __xdrWriter_h_

#include <rpc/types.h>
#include <rpc/xdr.h>

// TODO REMOVE THIS INCLUDE - TEMP HACK
#include "config.h"
#include "xdrWriter.h"

class XdrWriter
{
  private:
//    XDR *myXdr;
    FILE *myFile;

  public:
    XdrWriter(char* fileName, char* failureString, char* successString);
    ~XdrWriter();

    void static xdrWritePixel (ColPixel *col_pixel_p, XDR* myXdr, void (*ColourPalette) (float value, float col[]));
};

#endif //__xdrWriter_h_