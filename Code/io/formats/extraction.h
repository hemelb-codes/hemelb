#ifndef HEMELB_IO_FORMATS_EXTRACTION_H
#define HEMELB_IO_FORMATS_EXTRACTION_H

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      class Extraction
      {
        public:
          /**
           * Magic number to identify extraction data files.
           * ASCII for 'xtrx'
           */
          enum
          {
            MagicNumber = 0x78747278
          };

          /**
           * The version number of the file format.
           */
          enum
          {
            VersionNumber = 1
          };
      };
    }
  }
}

#endif /* HEMELB_IO_FORMATS_EXTRACTION_H */
