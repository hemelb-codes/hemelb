// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_MPIINCLUDE_H
#define HEMELB_MPIINCLUDE_H

#include <mpi.h>

namespace hemelb
{
  /*
   * There are templated functions for getting the MPI_Datatype
   * corresponding to a type or variable below. Use as:
   *     Foo bar;
   *     MPI_Datatype tp = MpiDataType(bar);
   *     // or alternatively
   *     tp = MpiDataType<Foo>();
   *
   * The generic version, for custom classes, is implemented using traits.
   * You must specialize the template, e.g.
   *
   *     template<>
   *     MPI_Datatype hemelb::MpiDataTypeTraits<Foo>::RegisterMpiDataType()
   *     {
   *       // Create the type
   *       MPI_Type_commit(&type);
   *       return type;
   *     }
   *
   * Built-in MPI types have specializations defined in the implementation.
   */

  template<typename T>
  class MpiDataTypeTraits
  {
    public:
      inline static const MPI_Datatype& GetMpiDataType()
      {
        if (mpiType == MPI_DATATYPE_NULL)
        {
          mpiType = RegisterMpiDataType();
        }
        return mpiType;
      }
    protected:
      static MPI_Datatype mpiType;
      static MPI_Datatype RegisterMpiDataType();
      template <class U> friend class MpiDataTypeTraits;
  };

  // Define the initial value
  template<typename T>
  MPI_Datatype MpiDataTypeTraits<T>::mpiType = MPI_DATATYPE_NULL;

  // Generic getters of data type info
  template<typename T>
  inline const MPI_Datatype& MpiDataType()
  {
    return MpiDataTypeTraits<T>::GetMpiDataType();
  }
  template<typename T>
  inline const MPI_Datatype& MpiDataType(const T&)
  {
    return MpiDataTypeTraits<T>::GetMpiDataType();
  }

}
#endif // HEMELB_MPIINCLUDE_H
