// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <fstream>
#include <memory>

#include "io/writers/AsciiFileWriter.h"

namespace hemelb::io
{

        // Implement a constructor that opens the file and creates the Xdr
        // object to write to it.
        AsciiFileWriter::AsciiFileWriter(std::string const& fileName)
        {
          outStream = std::make_unique<std::ofstream>(fileName);
        }

}
