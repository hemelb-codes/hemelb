#ifndef HEMELB_IO_WRITERS_DOUBLEOFFSETWRITER_H
#define HEMELB_IO_WRITERS_DOUBLEOFFSETWRITER_H

#include "io/writers/Writer.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      template<typename OutputType>
      class DoubleOffsetWriter
      {
        public:
          DoubleOffsetWriter(double offset, Writer& internalWriter) :
              offset(offset), internalWriter(internalWriter)
          {

          }

          DoubleOffsetWriter& operator<<(double const & value)
          {
            internalWriter << static_cast<OutputType>(value - offset);
            return *this;
          }

        private:
          // The offset to use when writing.
          double offset;
          // The actual writer object to pass through to.
          Writer& internalWriter;

      };
    }
  }
}

#endif /* HEMELB_IO_WRITERS_DOUBLEOFFSETWRITER_H */
