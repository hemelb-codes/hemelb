// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
