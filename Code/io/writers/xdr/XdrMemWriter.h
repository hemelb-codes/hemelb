// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRMEMWRITER_H
#define HEMELB_IO_WRITERS_XDR_XDRMEMWRITER_H

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
	// Simple wrapper to allow construction from pointer-length
	// arguments.
	class XdrMemWriter : public XdrMetaWriter<char*> {
	public:
	  XdrMemWriter(char* dataBuffer, unsigned int dataLength);
	};
      }
    }
  }
}
#endif // HEMELB_IO_WRITERS_XDR_XDRMEMWRITER_H
