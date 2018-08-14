
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

/*
 * XdrMemReader.cc
 *
 *  Created on: Oct 25, 2010
 *      Author: rupert
 */

#include "io/writers/xdr/XdrMemReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
        // Constructor to create an Xdr object based on a memory buffer
        XdrMemReader::XdrMemReader(const char* dataBuffer, unsigned int dataLength)
        {
          xdrmem_create(&mXdr, const_cast<char*>(dataBuffer), dataLength, XDR_DECODE);
        }

	XdrMemReader::XdrMemReader(const std::vector<char>& dataVec)
	{
	  xdrmem_create(&mXdr, const_cast<char*>(dataVec.data()), dataVec.size(), XDR_DECODE);
	}

      } // namespace xdr
    } // namespace writers
  }
}
