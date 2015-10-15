//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_IO_HDF5_TYPES_H
#define HEMELB_IO_HDF5_TYPES_H

#include <hdf5.h>

namespace hemelb
{
  namespace io
  {
    namespace hdf5
    {
      template<class T>
      struct Types
      {
      };

      template<>
      struct Types<char>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<char>::NATIVE = H5T_NATIVE_CHAR;

      template<>
      struct Types<double>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<double>::NATIVE = H5T_NATIVE_DOUBLE;

      template<>
      struct Types<float>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<float>::NATIVE = H5T_NATIVE_FLOAT;

      template<>
      struct Types<bool>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<bool>::NATIVE = H5T_NATIVE_HBOOL;

      template<>
      struct Types<long double>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<long double>::NATIVE = H5T_NATIVE_LDOUBLE;

      template<>
      struct Types<long>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<long>::NATIVE = H5T_NATIVE_LONG;

      template<>
      struct Types<long long>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<long long>::NATIVE = H5T_NATIVE_LLONG;

      template<>
      struct Types<short>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<short>::NATIVE = H5T_NATIVE_SHORT;

      template<>
      struct Types<unsigned long>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<unsigned long>::NATIVE = H5T_NATIVE_ULONG;

      template<>
      struct Types<unsigned long long>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<unsigned long long>::NATIVE = H5T_NATIVE_ULLONG;

      template<>
      struct Types<signed char>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<signed char>::NATIVE = H5T_NATIVE_SCHAR;

      template<>
      struct Types<unsigned char>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<unsigned char>::NATIVE = H5T_NATIVE_UCHAR;

      template<>
      struct Types<unsigned short>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<unsigned short>::NATIVE = H5T_NATIVE_USHORT;

      template<>
      struct Types<int>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<int>::NATIVE = H5T_NATIVE_INT;

      template<>
      struct Types<unsigned int>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<unsigned int>::NATIVE = H5T_NATIVE_UINT;

      template<>
      struct Types<char *>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<char *>::NATIVE = H5T_C_S1;

      template<>
      struct Types<std::string>
      {
          static const hid_t NATIVE;
      };
      const hid_t Types<std::string>::NATIVE = H5T_C_S1;
    }
  }
}

#endif  // HEMELB_IO_HDF5_TYPES_H
