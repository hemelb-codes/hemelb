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
  std::ostream * outStream;
  
  // Template for the simple write methods.
  template <typename T>
  void _write(T const & value) {
    *outStream << value;
  }
  
  // These necessary since can't override a virtual method with a
  // template member.
  void _write(int const & value);
  void _write(double const & value);
  void _write(float const & value);
  void _write(short const & value);
  void _write(unsigned int const & value);

};

#endif //__asciiStreamWriter_h_
