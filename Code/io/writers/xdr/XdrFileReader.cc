/*
 * XdrFileReader.cc
 *
 *  Created on: Oct 25, 2010
 *      Author: rupert
 */

#include "io/writers/xdr/XdrFileReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Constructor to create an Xdr object based on a file.
        XdrFileReader::XdrFileReader(FILE* xdrFile)
        {
          xdrstdio_create(&mXdr, xdrFile, XDR_DECODE);
        }

      } // namespace name

    } // namespace writers

  }

}
