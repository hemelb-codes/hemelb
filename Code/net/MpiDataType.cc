
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/MpiDataType.h"

namespace hemelb
{
  namespace net
  {
    // Specializations of the above getters for built in types.
    // These mappings are taken directly from the MPI standard version 2.2,
    // table 3.2 on page 28.
    template<>
    MPI_Datatype MpiDataTypeTraits<char>::RegisterMpiDataType()
    {
      return MPI_CHAR;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<int16_t>::RegisterMpiDataType()
    {
      return MPI_SHORT;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<int32_t>::RegisterMpiDataType()
    {
      return MPI_INT;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<int64_t>::RegisterMpiDataType()
    {
      return MPI_LONG_LONG;
    }
#ifdef __APPLE__
    template<>
    MPI_Datatype MpiDataTypeTraits<size_t>::RegisterMpiDataType()
    {
      return MPI_UNSIGNED_LONG;
    }
#endif
    template<>
    MPI_Datatype MpiDataTypeTraits<signed char>::RegisterMpiDataType()
    {
      return MPI_SIGNED_CHAR;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<unsigned char>::RegisterMpiDataType()
    {
      return MPI_UNSIGNED_CHAR;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<uint16_t>::RegisterMpiDataType()
    {
      return MPI_UNSIGNED_SHORT;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<uint32_t>::RegisterMpiDataType()
    {
      return MPI_UNSIGNED;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<uint64_t>::RegisterMpiDataType()
    {
      return MPI_UNSIGNED_LONG_LONG;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<float>::RegisterMpiDataType()
    {
      return MPI_FLOAT;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<double>::RegisterMpiDataType()
    {
      return MPI_DOUBLE;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<long double>::RegisterMpiDataType()
    {
      return MPI_LONG_DOUBLE;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<wchar_t>::RegisterMpiDataType()
    {
      return MPI_WCHAR;
    }
  }
}
