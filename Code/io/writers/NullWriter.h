// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_NULLWRITER_H
#define HEMELB_IO_WRITERS_NULLWRITER_H

#include <ostream>

#include "io/writers/Writer.h"

namespace hemelb::io
{
        class NullWriter : public Writer
        {
          public:
            // Method to get the current position of writing in the stream.
            unsigned int getCurrentStreamPosition() const override;

            // Methods for formatting control
            void writeFieldSeparator() override;
            void writeRecordSeparator() override;

          protected:

            template<typename T>
            void _write(T const & value)
            {
              // PASS, DO NOTHING
            }

            // These necessary since can't override a virtual method with a
            // template member.
            void _write(std::int16_t const& intToWrite) override;
            void _write(std::uint16_t const& uIntToWrite) override;
            void _write(std::int32_t const& intToWrite) override;
            void _write(std::uint32_t const& uIntToWrite) override;
            void _write(std::int64_t const& intToWrite) override;
            void _write(std::uint64_t const& uIntToWrite) override;

            void _write(double const& doubleToWrite) override;
            void _write(float const& floatToWrite) override;

            void _write(const std::string& floatToWrite) override;
        };

}
#endif // HEMELB_IO_WRITERS_NULLWRITER_H
