
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <fstream>

#include "io/writers/ascii/AsciiFileWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace ascii
      {

        // Implement a constructor that opens the file and creates the Xdr
        // object to write to it.
        AsciiFileWriter::AsciiFileWriter(std::string fileName)
        {
          outStream = new std::ofstream(fileName.c_str());
        }

        // Destructor closes the file and cleans up the stream object.
        AsciiFileWriter::~AsciiFileWriter()
        {
          if (outStream != NULL)
          {
            // Since we know the member std::ostream* outStream is actually an
            // instance of std::ofstream, cast it to that so we can close it.
            std::ofstream* outFileStream = static_cast<std::ofstream *>(outStream);

            outFileStream->close();
            delete outFileStream;
          }
        }
      } // namespace ascii
    } // namespace writers
  }
}
