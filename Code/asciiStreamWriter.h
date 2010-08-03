#ifndef __asciiStreamWriter_h_
#define __asciiStreamWriter_h_

#include <ostream>

#include "writer.h"

class AsciiStreamWriter : public Writer {
 public:
  // Method to get the current position of writing in the stream.
  unsigned int getCurrentStreamPosition() const;
  
  // Methods for formatting control
  void writeFieldSeparator();
  void writeRecordSeparator();

 protected:
  // Constructor and member variable. This class should never be
  // instantiated directly, so the constructor is only available to
  // subclasses.
  std::ostream * outStream;
  //AsciiStreamWriter();
  
  // Template for the simple write methods.
  template <typename T> void _write(T const & value);

};

#endif //__asciiStreamWriter_h_
