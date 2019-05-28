// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/writers/xdr/XdrFileWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

	// Use base class special constructor
	// Functor 1 opens the file
	// Functor 2 makes an iterator that writes to it
        XdrFileWriter::XdrFileWriter(const std::string& fileName) :
	  base([&fileName]() {
		 return std::ofstream(fileName);
	       },
	       [](std::ofstream& fb){
		 return std::ostreambuf_iterator<char>(fb);
	       })
	{
	}

      }
    }
  }
}
