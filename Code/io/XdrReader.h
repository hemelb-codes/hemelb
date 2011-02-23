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
        void readDouble(double& outDouble);
        void readFloat(float& outDouble);
        void readInt(int& outInt);
        void readUnsignedInt(unsigned int& outUInt);

        // Get the position in the stream.
        unsigned int GetPosition();

      protected:
        XdrReader();
        XDR mXdr;

    };
  }
}

#endif // HEMELB_IO_XDRREADER_H
