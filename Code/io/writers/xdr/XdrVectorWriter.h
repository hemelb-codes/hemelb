// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRVECTORWRITER_H
#define HEMELB_IO_WRITERS_XDR_XDRVECTORWRITER_H

#include "io/writers/xdr/XdrWriter.h"
#include <vector>

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

	// XDR encoder that will put its buffer in a vector so you
	// don't have to worry about managing memory.
	class XdrVectorWriter : public XdrWriter {
	public:
	  using VectorT = std::vector<char>;

	  XdrVectorWriter();
	  virtual ~XdrVectorWriter();

	  const VectorT& GetBuf() const;

	private:
	  template<typename, typename>
	  friend struct XdrVectorWriterHelper;
	  VectorT buffer;
	  FILE* file;
	};

      }
    }
  }
}

#endif
