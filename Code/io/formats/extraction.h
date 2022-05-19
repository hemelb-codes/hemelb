// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_EXTRACTION_H
#define HEMELB_IO_FORMATS_EXTRACTION_H

namespace hemelb::io::formats::extraction
{
  // Magic number to identify extraction data files.
  // ASCII for 'xtr' + EOF
  enum {
    MagicNumber = 0x78747204
  };

  // The version number of the file format.
  enum {
    VersionNumber = 5
  };

  // The length of the main header. Made up of:
  // uint - HemeLbMagicNumber
  // uint - ExtractionMagicNumber
  // uint - Format version number
  // double - Voxel size (metres)
  // double x 3 - Origin x,y,z components (metres)
  // uhyper - Total number of sites
  // uint - Field count
  // uint - Length of the field header that follows
  enum {
    MainHeaderLength = 60
    //!< MainHeaderLength
  };

  // Each field header is:
  // string - name
  // uint32 - number of elements
  // uint32 - type code
  // uint32 - number of offsets (valid values are {0, 1, n elem})
  // type[n offsets] - offsets (n offsets items of type implied above)
  enum class TypeCode : std::uint32_t {
    FLOAT,
    DOUBLE,
    INT32,
    UINT32,
    INT64,
    UINT64,
  };

  // Compute the length of data written by XDR for a given string.
  inline size_t GetStoredLengthOfString(std::string const& str)
  {
    // XDR pads up to the nearest multiple of four bytes and also
    // stores the length of the string.
    size_t len = str.length();
    size_t mod = len % 4;
    if (mod) {
      len += 4 - mod;
    }
    return len + 4;
  }

  // Compute the length (in bytes) of a field header
  inline size_t GetFieldHeaderLength(std::string const& name, std::uint32_t noff, TypeCode tc) {
    size_t len = GetStoredLengthOfString(name);
    len += 4;  // number of elements
    len += 4;  // type code
    len += 4;  // number of offsets
    size_t elemsize = (tc == TypeCode::FLOAT || tc == TypeCode::INT32 || tc == TypeCode::UINT32) ? 4 : 8;
    len += elemsize * noff; // offsets
    return len;
  }

}
#endif // HEMELB_IO_FORMATS_EXTRACTION_H
