// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_ASCIIFILEWRITER_H
#define HEMELB_IO_WRITERS_ASCIIFILEWRITER_H

#include <iostream>

#include "io/writers/AsciiStreamWriter.h"

namespace hemelb::io
{
        // Class to write a file. The actual write functions are implemented
        // in the base class.
        class AsciiFileWriter : public AsciiStreamWriter
        {

          public:
            explicit AsciiFileWriter(std::string const& fileName);
            ~AsciiFileWriter() override = default;

        };
}
#endif // HEMELB_IO_WRITERS_ASCIIFILEWRITER_H
