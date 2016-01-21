
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_EXTRACTION_H
#define HEMELB_IO_FORMATS_EXTRACTION_H

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      namespace extraction
      {
        /**
         * Magic number to identify extraction data files.
         * ASCII for 'xtr' + EOF
         */
        enum
        {
          MagicNumber = 0x78747204
        };

        /**
         * The version number of the file format.
         */
        enum
        {
          VersionNumber = 4
        };

        /**
         * The length of the main header. Made up of:
         * uint - HemeLbMagicNumber
         * uint - ExtractionMagicNumber
         * uint - Format version number
         * double - Voxel size (metres)
         * double x 3 - Origin x,y,z components (metres)
         * uhyper - Total number of sites
         * uint - Field count
         * uint - Length of the field header that follows
         */
        enum
        {
          MainHeaderLength = 60
        //!< MainHeaderLength
        };

        /**
         * Compute the length of data written by XDR for a given string.
         * @param str
         * @return
         */
        inline size_t GetStoredLengthOfString(std::string str)
        {
          // XDR pads up to the nearest multiple of four bytes and also
          // stores the length of the string.
          size_t len = str.length();
          size_t mod = len % 4;
          if (mod)
          {
            len += 4 - mod;
          }
          return len + 4;
        }
      }
    }
  }
}
#endif /* HEMELB_IO_FORMATS_EXTRACTION_H */
