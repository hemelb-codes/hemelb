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
#include "Exception.h"

namespace hemelb
{
  namespace net
  {
    template<typename vectorType>
    MPI_Datatype GenerateTypeForVector()
    {
      const int typeCount = 1;
      int blocklengths[typeCount] = { 3 };

      MPI_Datatype types[typeCount] = { net::MpiDataType<vectorType>() };

      MPI_Aint displacements[typeCount] = { 0 };

      MPI_Datatype ret;

      HEMELB_MPI_CALL(MPI_Type_create_struct, (typeCount, blocklengths, displacements, types, &ret));

      HEMELB_MPI_CALL(MPI_Type_commit, (&ret));
      return ret;
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<float> >::RegisterMpiDataType()
    {
      return GenerateTypeForVector<float>();
    }

    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<site_t> >::RegisterMpiDataType()
    {
      return GenerateTypeForVector<site_t>();
    }

    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::util::Vector3D<distribn_t> >::RegisterMpiDataType()
    {
      return GenerateTypeForVector<distribn_t>();
    }
  }
  namespace util
  {
    namespace
    {
      void DefaultHandlerFunction(int direction)
      {
        // TODO need to find a way of handling this case better.
        throw Exception() << "Failed while accessing a direction in Vector3D.";
      }
    }
    Vector3DBase::HandlerFunction* Vector3DBase::handler = DefaultHandlerFunction;

  }
}
