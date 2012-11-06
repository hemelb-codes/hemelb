// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "io/writers/xdr/XdrReader.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        XdrReader::XdrReader()
        {
        }

        // Functions to read out the next bit of the file as a certain type.
        bool XdrReader::readDouble(double& outDouble)
        {
          return xdr_double(&mXdr, &outDouble);
        }

        bool XdrReader::readFloat(float& outDouble)
        {
          return xdr_float(&mXdr, &outDouble);
        }

        bool XdrReader::readInt(int& outInt)
        {
          return xdr_int(&mXdr, &outInt);
        }

        bool XdrReader::readUnsignedInt(unsigned int& outUInt)
        {
          return xdr_u_int(&mXdr, &outUInt);
        }

        bool XdrReader::readUnsignedLong(uint64_t& outULong)
        {
          u_quad_t temporary;
          bool ret = xdr_uint64_t(&mXdr, &temporary);
          outULong = temporary;
          return ret;
        }

        unsigned int XdrReader::GetPosition()
        {
          return xdr_getpos(&mXdr);
        }

        // Returns false on failure
        bool XdrReader::SetPosition(unsigned int iPosition)
        {
          return xdr_setpos(&mXdr, iPosition);
        }

        // Destructor to get rid of any resources used by the Xdr object. This class doesn't create
        // the file object, so it doesn't free it either.
        XdrReader::~XdrReader()
        {
          xdr_destroy(&mXdr);
        }

      } // namespace xdr
    } // namespace writers
  }
}
