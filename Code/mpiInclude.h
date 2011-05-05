#ifndef HEMELB_MPIINCLUDE_H
#define HEMELB_MPIINCLUDE_H

#ifdef XT3
#include <mpi.h>
#else
#include "mpi.h"
#endif

/*
 * Templated function for getting the MPI_Datatype corresponding to a type.
 * The generic version, for custom classes, expects that the class T has a
 * static member function const MPI_DataType T::GetMpiDataType().
 *
 * Builtin MPI types have specialisations defined.
 */
template<typename T>
inline const MPI_Datatype MpiDataType()
{
  return T::GetMpiDataType();
}

// These mappings are taken directly from the MPI standard version 2.2,
// table 3.2 on page 28.
template<>
inline const MPI_Datatype MpiDataType<char> ()
{
  return MPI_CHAR;
}
template<>
inline const MPI_Datatype MpiDataType<signed short int> ()
{
  return MPI_SHORT;
}
template<>
inline const MPI_Datatype MpiDataType<signed int> ()
{
  return MPI_INT;
}
template<>
inline const MPI_Datatype MpiDataType<signed long int> ()
{
  return MPI_LONG;
}
// Strictly, C++ doesn't have long long
//template<>
//inline const MPI_Datatype MpiDataType<signed long long int> ()
//{
//  return MPI_LONG_LONG_INT;
//}
template<>
inline const MPI_Datatype MpiDataType<signed char> ()
{
  return MPI_SIGNED_CHAR;
}
template<>
inline const MPI_Datatype MpiDataType<unsigned char> ()
{
  return MPI_UNSIGNED_CHAR;
}
template<>
inline const MPI_Datatype MpiDataType<unsigned short int> ()
{
  return MPI_UNSIGNED_SHORT;
}
template<>
inline const MPI_Datatype MpiDataType<unsigned int> ()
{
  return MPI_UNSIGNED;
}
template<>
inline const MPI_Datatype MpiDataType<unsigned long int> ()
{
  return MPI_UNSIGNED_LONG;
}
// Strictly, C++ doesn't have long long
//template<>
//inline const MPI_Datatype MpiDataType<unsigned long long int> ()
//{
//  return MPI_UNSIGNED_LONG_LONG;
//}
template<>
inline const MPI_Datatype MpiDataType<float> ()
{
  return MPI_FLOAT;
}
template<>
inline const MPI_Datatype MpiDataType<double> ()
{
  return MPI_DOUBLE;
}
template<>
inline const MPI_Datatype MpiDataType<long double> ()
{
  return MPI_LONG_DOUBLE;
}
template<>
inline const MPI_Datatype MpiDataType<wchar_t> ()
{
  return MPI_WCHAR;
}


#endif // HEMELB_MPIINCLUDE_H
