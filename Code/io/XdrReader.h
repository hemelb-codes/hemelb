#ifndef HEMELB_IO_XDRREADER_H
#define HEMELB_IO_XDRREADER_H

#include <rpc/types.h>
#include <rpc/xdr.h>

namespace hemelb
{
  namespace io
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

        // Get the position in the stream.
        unsigned int GetPosition();
        bool SetPosition(unsigned int iPosition);

      protected:
        XdrReader();
        XDR mXdr;

    };
  }
}

#endif // HEMELB_IO_XDRREADER_H
