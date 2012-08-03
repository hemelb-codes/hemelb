// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

        // Class for writing to an Xdr object in memory. Derives from a
        // base XdrWriter which implements all the writing functions.

        class XdrMemWriter : public XdrWriter
        {
          public:
            // Constructor and destructor for the in-memory Xdr writer.
            XdrMemWriter(char* dataBuffer, unsigned int dataLength);
            ~XdrMemWriter();
        };

      } // namespace xdr
    } // namespace writers
  }
}
#endif // HEMELB_IO_WRITERS_XDR_XDRMEMWRITER_H
