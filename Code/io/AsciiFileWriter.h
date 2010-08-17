#ifndef __io_AsciiFileWriter_h_
#define __io_AsciiFileWriter_h_

#include <iostream>

#include "io/AsciiStreamWriter.h"

namespace io {
  // Class to write a file. The actual write functions are implemented
  // in the base class.
  class AsciiFileWriter : public AsciiStreamWriter {
  
  public:
    AsciiFileWriter(char* fileName);
    ~AsciiFileWriter();
  
  };
  
}

#endif //__io_AsciiFileWriter_h_
