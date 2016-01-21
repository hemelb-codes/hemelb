
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/writers/xdr/XdrMemWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Constructor for a Xdr writer held in a memory buffer.
        XdrMemWriter::XdrMemWriter(char* dataBuffer, unsigned int dataLength)
        {
          xdrmem_create(&mXdr, dataBuffer, dataLength, XDR_ENCODE);
        }

        // Destructor for the class.
        XdrMemWriter::~XdrMemWriter()
        {
          xdr_destroy(&mXdr);
        }

      } // namespace xdr
    } // namespace writers
  }
}
