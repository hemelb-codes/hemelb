#include "io/writers/Writer.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
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

    } // namespace writer
  }
}
