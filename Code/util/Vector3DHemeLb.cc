// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "util/Vector3D.h"
#include "net/mpi.h"
#include "units.h"
#include <stdexcept>
#include "log/Logger.h"

namespace hemelb
{
  template<typename vectorType>
  MPI_Datatype GenerateTypeForVector()
  {
    const int typeCount = 1;
    int blocklengths[typeCount] = { 3 };

    MPI_Datatype types[typeCount] = { net::MpiDataType<vectorType> () };

    MPI_Aint displacements[typeCount] = { 0 };

    MPI_Datatype ret;

    MPI_Type_struct(typeCount, blocklengths, displacements, types, &ret);

    MPI_Type_commit(&ret);
    return ret;
  }

  template<>
  MPI_Datatype net::MpiDataTypeTraits<hemelb::util::Vector3D<float> >::RegisterMpiDataType()
  {
    return GenerateTypeForVector<float> ();
  }

  template<>
  MPI_Datatype net::MpiDataTypeTraits<hemelb::util::Vector3D<site_t> >::RegisterMpiDataType()
  {
    return GenerateTypeForVector<site_t> ();
  }

  template<>
   MPI_Datatype net::MpiDataTypeTraits<hemelb::util::Vector3D<distribn_t> >::RegisterMpiDataType()
   {
     return GenerateTypeForVector<distribn_t> ();
   }
  namespace util
  {
    namespace
    {
      void DefaultHandlerFunction(int direction)
      {
        // TODO need to find a way of handling this case better.
        throw std::runtime_error("Failed while accessing a direction in Vector3D.");
      }
    }
    Vector3DBase::HandlerFunction* Vector3DBase::handler = DefaultHandlerFunction;

  }
}
