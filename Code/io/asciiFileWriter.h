#ifndef __io_asciiFileWriter_h_
#define __io_asciiFileWriter_h_

#include <iostream>

#include "io/asciiStreamWriter.h"

namespace io {
  // Class to write a file. The actual write functions are implemented
  // in the base class.
  class AsciiFileWriter : public AsciiStreamWriter {
  
  public:
    AsciiFileWriter(char* fileName);
    ~AsciiFileWriter();
  
  };
  
}

#endif //__io_asciiFileWriter_h_
