#ifndef HEMELB_IO_FORMATS_H
#define HEMELB_IO_FORMATS_H

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      // This is the generic HemeLB binary file identifier. It should be the
      // first 4 bytes of every (binary) file used for IO. It should be then
      // followed by another 4 bytes identifying the particular type/version,
      // that number being terminated by the EOF character (0x04)
      enum
      {
        // ASCII for hlb!
        HemeLbMagicNumber = 0x686c6221
      };
    }
  }

}
#endif // HEMELB_IO_FORMATS_H
