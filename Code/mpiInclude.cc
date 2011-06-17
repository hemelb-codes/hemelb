#include <stdint.h>

#include "mpiInclude.h"

namespace hemelb
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
  MPI_Datatype MpiDataTypeTraits<signed short int>::RegisterMpiDataType()
  {
    return MPI_SHORT;
  }
  template<>
  MPI_Datatype MpiDataTypeTraits<signed int>::RegisterMpiDataType()
  {
    return MPI_INT;
  }
  template<>
  MPI_Datatype MpiDataTypeTraits<signed long int>::RegisterMpiDataType()
  {
    return MPI_LONG_LONG;
  }
  // Strictly, C++ doesn't have long long
  //template<>
  //MPI_Datatype MpiDataTypeTraits<int64_t>::RegisterMpiDataType()
  //{
  //  return MPI_LONG_LONG;
  //}
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
  MPI_Datatype MpiDataTypeTraits<unsigned short int>::RegisterMpiDataType()
  {
    return MPI_UNSIGNED_SHORT;
  }
  template<>
  MPI_Datatype MpiDataTypeTraits<unsigned int>::RegisterMpiDataType()
  {
    return MPI_UNSIGNED;
  }
  template<>
  MPI_Datatype MpiDataTypeTraits<unsigned long int>::RegisterMpiDataType()
  {
    return MPI_UNSIGNED_LONG_LONG;
  }
  // Strictly, C++ doesn't have long long
  //template<>
  //MPI_Datatype MpiDataTypeTraits<unsigned long long int>::RegisterMpiDataType()
  //{
  //  return MPI_UNSIGNED_LONG_LONG;
  //}
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
