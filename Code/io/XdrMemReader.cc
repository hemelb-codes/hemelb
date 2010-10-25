/*
 * XdrFileReader.cc
 *
 *  Created on: Oct 25, 2010
 *      Author: rupert
 */

#include "XdrMemReader.h"

namespace hemelb
{

  namespace io
  {
    // Constructor to create an Xdr object based on a memory buffer
    XdrMemReader::XdrMemReader(char* dataBuffer, unsigned int dataLength)
    {
      mXdr = new XDR;
      xdrmem_create(mXdr, dataBuffer, dataLength, XDR_DECODE);
    }


  }

}
