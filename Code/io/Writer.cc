#include "io/Writer.h"
#include "vis/ColPixel.h"

namespace hemelb
{
  namespace io
  {

    Writer::Writer()
    {
    }

    Writer::~Writer()
    {
      // Pure virtual destructor
    }

    Writer& Writer::operator<<(enum Writer::Separator const & value)
    {
      writeRecordSeparator();
      return *this;
    }

  }
}
