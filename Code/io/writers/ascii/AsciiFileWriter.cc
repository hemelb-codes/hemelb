// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
          if (outStream != nullptr)
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
