#ifndef HEMELB_IO_ASCIISTREAMWRITER_H
#define HEMELB_IO_ASCIISTREAMWRITER_H

#include <ostream>

#include "io/Writer.h"

namespace hemelb
{
  namespace io
  {

    class AsciiStreamWriter : public Writer
    {
      public:
        // Method to get the current position of writing in the stream.
        unsigned int getCurrentStreamPosition() const;

        // Methods for formatting control
        void writeFieldSeparator();
        void writeRecordSeparator();

      protected:
        std::ostream * outStream;

        template<typename T>
        void _write(T const & value)
        {
          *outStream << value;
        }

        // These necessary since can't override a virtual method with a
        // template member.
        void _write(int16_t const& intToWrite);
        void _write(uint16_t const& uIntToWrite);
        void _write(int32_t const& intToWrite);
        void _write(uint32_t const& uIntToWrite);
        void _write(int64_t const& intToWrite);
        void _write(uint64_t const& uIntToWrite);

        void _write(double const& doubleToWrite);
        void _write(float const& floatToWrite);

    };
  }
}
#endif // HEMELB_IO_ASCIISTREAMWRITER_H
