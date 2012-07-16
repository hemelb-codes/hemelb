// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_WRITERS_ASCII_ASCIIFILEWRITER_H
#define HEMELB_IO_WRITERS_ASCII_ASCIIFILEWRITER_H

#include <iostream>

#include "io/writers/ascii/AsciiStreamWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace ascii
      {
        // Class to write a file. The actual write functions are implemented
        // in the base class.
        class AsciiFileWriter : public AsciiStreamWriter
        {

          public:
            AsciiFileWriter(std::string fileName);
            ~AsciiFileWriter();

        };
      } // namespace ascii
    } // namespace writers
  }
}
#endif // HEMELB_IO_WRITERS_ASCII_ASCIIFILEWRITER_H
