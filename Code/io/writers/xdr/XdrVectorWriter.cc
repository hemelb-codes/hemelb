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
	// The plan here it to use funopen (BSDs, ie macOS) or
	// fcookieopen (Linux with glibc) to create a FILE* object
	// backed by a std::vector

	// This class actually implements the writing to the
	// std::vector. It has to be templated on the different
	// signatures required by funopen and cookie_io_functions_t
	template <typename nwritten_t, typename nbyte_t>
	struct XdrVectorWriterHelper {
	  // The equivalent of fwrite(FILE*, ...) but we know that the
	  // "cookie" was the pointer to the writer passed in by
	  // fopen_wrapper.
	  static nwritten_t fwrite(void* cookie, const char* data, nbyte_t nbyte) {
	    auto owner = reinterpret_cast<XdrVectorWriter*>(cookie);
	    auto& buf = owner->buffer;
	    buf.reserve(buf.size() + nbyte);
	    buf.insert(buf.end(), data, data + nbyte);
	    return nbyte;
	  }
	};

	// Make the relevant syscall to create our fake FILE*
	static FILE* fopen_wrapper(XdrVectorWriter* owner) {
#if defined(HAVE_FUNOPEN)
	  // macOS and other BSDs
	  return funopen(owner,
			 nullptr,
			 XdrVectorWriterHelper<int, int>::fwrite,
			 nullptr,
			 nullptr);
#elif defined(HAVE_FOPENCOOKIE)
	  // Linux with GNU standard library
	    const char* mode = "w";
	    cookie_io_functions_t fns = {nullptr,
					 XdrVectorWriterHelper<ssize_t, size_t>::fwrite,
					 nullptr,
					 nullptr};
	    return fopencookie(owner, mode, fns);
#else
#error "Don't know how to create a file-like interface on this platform"
	    return nullptr;
#endif
	}
	
	XdrVectorWriter::XdrVectorWriter() : buffer(), file(nullptr) {
	  file = fopen_wrapper(this);
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
