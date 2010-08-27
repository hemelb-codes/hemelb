#ifndef HEMELB_IO_ASCIIFILEWRITER_H
#define HEMELB_IO_ASCIIFILEWRITER_H

#include <iostream>

#include "io/AsciiStreamWriter.h"

namespace hemelb
{
  namespace io
  {
    
    // Class to write a file. The actual write functions are implemented
    // in the base class.
    class AsciiFileWriter : public AsciiStreamWriter {
  
    public:
      AsciiFileWriter(char* fileName);
      ~AsciiFileWriter();
  
    };
  
  }

}
#endif // HEMELB_IO_ASCIIFILEWRITER_H
