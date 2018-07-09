
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDR_XDRREADER_H
#define HEMELB_IO_WRITERS_XDR_XDRREADER_H

#if HEMELB_HAVE_CSTDINT
# include <cstdint>
#else
# include <stdint.h>
#endif
#include <rpc/types.h>
#include <rpc/xdr.h>

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {
        // Class to read an Xdr-style file from disk.
        class XdrReader
        {
          public:
            // destructor.
            virtual ~XdrReader();

            // Functions for reading the next bit of the stream.
            bool readDouble(double& outDouble);
            bool readFloat(float& outDouble);
            bool readInt(int& outInt);
            bool readUnsignedInt(unsigned int& outUInt);
            bool readUnsignedLong(uint64_t& outULong);

            // Get the position in the stream.
            unsigned int GetPosition();
            bool SetPosition(unsigned int iPosition);

          protected:
            XdrReader();
            XDR mXdr;

        };

      } // namespace xdr
    } // namespace writers
  }
}

#endif // HEMELB_IO_WRITERS_XDR_XDRREADER_H
