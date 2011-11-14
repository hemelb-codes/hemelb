/*
 * XdrFileReader.cc
 *
 *  Created on: Oct 25, 2010
 *      Author: rupert
 */

#include "XdrFileReader.h"

namespace hemelb
{

  namespace io
  {
    // Constructor to create an Xdr object based on a file.
    XdrFileReader::XdrFileReader(FILE* xdrFile)
    {
      xdrstdio_create(&mXdr, xdrFile, XDR_DECODE);
    }

  }

}
