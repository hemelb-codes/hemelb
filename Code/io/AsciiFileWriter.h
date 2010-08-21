#ifndef HEME_IO_ASCIIFILEWRITER_H
#define HEME_IO_ASCIIFILEWRITER_H

#include <iostream>

#include "io/AsciiStreamWriter.h"

namespace heme {
  namespace io {
    
    // Class to write a file. The actual write functions are implemented
    // in the base class.
    class AsciiFileWriter : public AsciiStreamWriter {
  
    public:
      AsciiFileWriter(char* fileName);
      ~AsciiFileWriter();
  
    };
  
  }

}
#endif // HEME_IO_ASCIIFILEWRITER_H
