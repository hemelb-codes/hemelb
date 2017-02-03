//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MPIDATATYPE_H
#define HEMELB_NET_MPIDATATYPE_H

#include <mpi.h>

#define HEMELB_MPI_TYPE_BEGIN(outType, Type, n) \
  MPI_Datatype outType = MPI_DATATYPE_NULL; \
  { \
  Type typeInstances[2]; \
  MPI_Aint instanceAddr; \
  HEMELB_MPI_CALL(MPI_Get_address, (typeInstances, &instanceAddr)); \
  const unsigned elementCount = n; \
  unsigned elementCounter = 0; \
  int elementBlockLengths[elementCount]; \
  MPI_Aint elementDisplacements[elementCount]; \
  MPI_Datatype elementTypes[elementCount]

#define HEMELB_MPI_TYPE_ADD_MEMBER_N(name, count) \
  if (elementCounter >= elementCount) throw ::hemelb::Exception() \
    << "Attempting to define more members than specified"; \
  elementBlockLengths[elementCounter] = count; \
  HEMELB_MPI_CALL(MPI_Get_address, \
		  (&(typeInstances->name), elementDisplacements + elementCounter) \
		  ); \
  elementDisplacements[elementCounter] -= instanceAddr; \
  elementTypes[elementCounter] = ::hemelb::net::MpiDataType(typeInstances->name); \
  ++elementCounter

#define HEMELB_MPI_TYPE_ADD_MEMBER(name) HEMELB_MPI_TYPE_ADD_MEMBER_N(name, 1)

#define HEMELB_MPI_TYPE_END(outType, Type) \
  if (elementCounter != elementCount) throw ::hemelb::Exception() \
    << "Error in type definition: only " << elementCounter << " elements added but specified " << elementCount; \
  MPI_Datatype tmp_type; \
  HEMELB_MPI_CALL(MPI_Type_create_struct, \
                  (elementCount, elementBlockLengths, elementDisplacements, elementTypes, &tmp_type) \
                  ); \
  MPI_Aint instanceSize; \
  HEMELB_MPI_CALL(MPI_Get_address, (&typeInstances[1], &instanceSize)); \
  instanceSize -= instanceAddr; \
  HEMELB_MPI_CALL(MPI_Type_create_resized, \
                  (tmp_type, 0, instanceSize, &outType) \
                  ); \
  HEMELB_MPI_CALL(MPI_Type_free, \
                  (&tmp_type) \
                  ); \
  }

namespace hemelb
{
  namespace net
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
     *     MPI_Datatype hemelb::net::MpiDataTypeTraits<Foo>::RegisterMpiDataType()
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
      private:
        static MPI_Datatype mpiType;
        static MPI_Datatype RegisterMpiDataType();
        template<class U> friend class MpiDataTypeTraits;
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
}
#endif // HEMELB_NET_MPIDATATYPE_H
