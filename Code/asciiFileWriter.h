#ifndef __asciiFileWriter_h_
#define __asciiFileWriter_h_

#include <iostream>

#include "asciiStreamWriter.h"

// Class to write a file. The actual write functions are implemented
// in the base class.
class AsciiFileWriter : public AsciiStreamWriter {
  
 public:
  AsciiFileWriter(char* fileName);
  ~AsciiFileWriter();
  
};

#endif //__asciiFileWriter_h_
