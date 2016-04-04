
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>

#include "io/writers/null/NullWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace null
      {

        void NullWriter::_write(int16_t const & value)
        {
          this->_write<int16_t> (value);
        }
        void NullWriter::_write(uint16_t const & value)
        {
          this->_write<uint16_t> (value);
        }
        void NullWriter::_write(int32_t const & value)
        {
          this->_write<int32_t> (value);
        }
        void NullWriter::_write(uint32_t const & value)
        {
          this->_write<uint32_t> (value);
        }
        void NullWriter::_write(int64_t const & value)
        {
          this->_write<int64_t> (value);
        }
        void NullWriter::_write(uint64_t const & value)
        {
          this->_write<uint64_t> (value);
        }
        void NullWriter::_write(double const & value)
        {
          this->_write<double> (value);
        }
        void NullWriter::_write(float const & value)
        {
          this->_write<float> (value);
        }

        void NullWriter::_write(const std::string& value)
        {
          this->_write<std::string> (value);
        }

      } // namespace ascii

    } // namespace writers
  }
}
