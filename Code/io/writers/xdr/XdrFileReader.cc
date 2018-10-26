
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

/*
 * XdrFileReader.cc
 *
 *  Created on: Oct 25, 2010
 *      Author: rupert
 */
#include <cassert>
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
        XdrFileReader::XdrFileReader(const std::string& fn) {
	  fh = std::fopen(fn.c_str(), "r");
	  assert(fh != nullptr);
          xdrstdio_create(&mXdr, fh, XDR_DECODE);
        }

	XdrFileReader::~XdrFileReader() {
	  std::fclose(fh);
	}
      } // namespace name

    } // namespace writers

  }

}
