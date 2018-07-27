// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <stdio.h>
#include "io/writers/xdr/XdrVectorWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

	// Helper struct with static methods to pass to funopen
	struct XdrVectorWriterHelper {
	  static FILE* fopen(XdrVectorWriter* owner) {
	    return funopen(owner, nullptr, XdrVectorWriterHelper::fwrite, nullptr, nullptr);
	  }

	  static int fwrite(void* cookie, const char* data, int nbyte) {
	    auto owner = reinterpret_cast<XdrVectorWriter*>(cookie);
	    auto& buf = owner->buffer;
	    buf.reserve(buf.size() + nbyte);
	    buf.insert(buf.end(), data, data + nbyte);
	    return nbyte;
	  }
	};

	XdrVectorWriter::XdrVectorWriter() : buffer(), file(nullptr) {
	  file = XdrVectorWriterHelper::fopen(this);
	  xdrstdio_create(&mXdr, file, XDR_ENCODE);
	}

	XdrVectorWriter::~XdrVectorWriter() {
	  xdr_destroy(&mXdr);
	  fclose(file);
	}

	auto XdrVectorWriter::GetBuf() const -> const VectorT& {
	  fflush(file);
	  return buffer;
	}

      } // namespace xdr
    }
  }
}
