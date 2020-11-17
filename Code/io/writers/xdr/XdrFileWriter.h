// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRFILEWRITER_H
#define HEMELB_IO_WRITERS_XDR_XDRFILEWRITER_H

#include <string>
#include <fstream>
#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

	using CharFileIter = std::ostreambuf_iterator<char>;

	// Wrapper that owns a filehandle and uses that as the sink for data
        class XdrFileWriter : public XdrMetaWriter<CharFileIter, std::ofstream>
	{
	  using base = XdrMetaWriter<CharFileIter, std::ofstream>;
  	public:
	  XdrFileWriter(const std::string& fileName);
	};

      }
    }
  }
}
#endif // HEMELB_IO_WRITERS_XDR_XDRFILEWRITER_H
