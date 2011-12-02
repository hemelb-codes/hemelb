#include <rpc/types.h>
#include <rpc/xdr.h>

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    // Functions to write simple types out to the Xdr stream.
    // Sadly templating is non-trivial due to XDR's naming.

    /* The const_cast in these is safe as the XDR writer won't
     * change the values but doesn't specify in the headers.
     */


    namespace writers
    {
      namespace xdr
      {
        /* The XDR implementation on BSD based machines uses
            * xdr_u_int??_t() to pack/unpack these types, while Linux
            * based machines have xdr_uint??_t(). Use the config macro to
            * create aliases on BSD.
            */
#ifdef HEMELB_CFG_ON_BSD
#define xdr_uint16_t xdr_u_int16_t
#define xdr_uint32_t xdr_u_int32_t
#define xdr_uint64_t xdr_u_int64_t
#endif //  HEMELB_CFG_ON_BSD
        void XdrWriter::_write(int16_t const& shortToWrite)
        {
          xdr_int16_t(&mXdr, const_cast<int16_t *>(&shortToWrite));
        }

        void XdrWriter::_write(uint16_t const& shortToWrite)
        {
          xdr_uint16_t(&mXdr, const_cast<uint16_t *>(&shortToWrite));
        }

        void XdrWriter::_write(int32_t const& value)
        {
          xdr_int(&mXdr, const_cast<int32_t *>(&value));
        }

        void XdrWriter::_write(uint32_t const& uIntToWrite)
        {
          xdr_uint32_t(&mXdr, const_cast<uint32_t *>(&uIntToWrite));
        }

        void XdrWriter::_write(int64_t const& longToWrite)
        {
          xdr_int64_t(&mXdr, const_cast<int64_t*>(&longToWrite));
        }

        void XdrWriter::_write(uint64_t const& longToWrite)
        {
          xdr_uint64_t(&mXdr, const_cast<uint64_t*>(&longToWrite));
        }

        void XdrWriter::_write(float const& floatToWrite)
        {
          xdr_float(&mXdr, const_cast<float *>(&floatToWrite));
        }

        void XdrWriter::_write(double const& doubleToWrite)
        {
          xdr_double(&mXdr, const_cast<double *>(&doubleToWrite));
        }

        // Method to get the current position in the stream.
        unsigned int XdrWriter::getCurrentStreamPosition() const
        {
          // The XDR function does not modify the instance, but is not const in the declaration.
          return xdr_getpos(const_cast<XDR*> (&mXdr));
        }

        // No field/record separators in XDR files
        void XdrWriter::writeFieldSeparator()
        {
        }
        void XdrWriter::writeRecordSeparator()
        {
        }

      } // namespace xdr
    } // namespace writers
  }
}
