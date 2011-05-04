#include <iostream>

#include "io/AsciiStreamWriter.h"

namespace hemelb
{
  namespace io
  {

    // Method to get the current position in the stream.
    unsigned int AsciiStreamWriter::getCurrentStreamPosition() const
    {
      return (unsigned int) outStream->tellp();
    }

    // No field/record separators in XDR files
    void AsciiStreamWriter::writeFieldSeparator()
    {
      *outStream << ' ';
    }
    void AsciiStreamWriter::writeRecordSeparator()
    {
      *outStream << std::endl;
    }

    void AsciiStreamWriter::_write(int16_t const & value)
    {
      this->_write<int16_t> (value);
    }
    void AsciiStreamWriter::_write(u_int16_t const & value)
    {
      this->_write<u_int16_t> (value);
    }
    void AsciiStreamWriter::_write(int32_t const & value)
    {
      this->_write<int32_t> (value);
    }
    void AsciiStreamWriter::_write(u_int32_t const & value)
    {
      this->_write<u_int32_t> (value);
    }
    void AsciiStreamWriter::_write(int64_t const & value)
    {
      this->_write<int64_t> (value);
    }
    void AsciiStreamWriter::_write(u_int64_t const & value)
    {
      this->_write<u_int64_t> (value);
    }
    void AsciiStreamWriter::_write(double const & value)
    {
      this->_write<double> (value);
    }
    void AsciiStreamWriter::_write(float const & value)
    {
      this->_write<float> (value);
    }
  }
}
