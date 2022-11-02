// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>

#include "io/writers/AsciiStreamWriter.h"

namespace hemelb::io
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
          this->_write<int16_t>(value);
        }
        void AsciiStreamWriter::_write(uint16_t const & value)
        {
          this->_write<uint16_t>(value);
        }
        void AsciiStreamWriter::_write(int32_t const & value)
        {
          this->_write<int32_t>(value);
        }
        void AsciiStreamWriter::_write(uint32_t const & value)
        {
          this->_write<uint32_t>(value);
        }
        void AsciiStreamWriter::_write(int64_t const & value)
        {
          this->_write<int64_t>(value);
        }
        void AsciiStreamWriter::_write(uint64_t const & value)
        {
          this->_write<uint64_t>(value);
        }
        void AsciiStreamWriter::_write(double const & value)
        {
          this->_write<double>(value);
        }
        void AsciiStreamWriter::_write(float const & value)
        {
          this->_write<float>(value);
        }

        void AsciiStreamWriter::_write(const std::string& value)
        {
          this->_write<std::string>(value);
        }

}
