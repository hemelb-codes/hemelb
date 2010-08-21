#ifndef HEME_IO_ASCIISTREAMWRITER_H
#define HEME_IO_ASCIISTREAMWRITER_H

#include <ostream>

#include "io/Writer.h"

namespace heme {
  namespace io {
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
  }
}
#endif // HEME_IO_ASCIISTREAMWRITER_H
