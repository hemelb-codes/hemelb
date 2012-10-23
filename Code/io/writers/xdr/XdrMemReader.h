// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRMEMREADER_H

#include <cstdio>

#include "io/writers/xdr/XdrReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        class XdrMemReader : public XdrReader
        {
          public:
            XdrMemReader(char* dataBuffer, unsigned int dataLength);

        };
      } // namespace xdr
    } // namespace writers
  }
}

#endif /* HEMELB_IO_XDRFILEREADER_H */
