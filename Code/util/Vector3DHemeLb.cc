// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "util/Vector3D.h"
#include "net/mpi.h"
#include "Exception.h"

namespace hemelb::net
{
    template<>
    MPI_Datatype MpiDataTypeTraits<util::Vector3D<float> >::RegisterMpiDataType()
    {
      return MpiDataTypeTraits<std::array<float, 3>>::RegisterMpiDataType();
    }

    template<>
    MPI_Datatype MpiDataTypeTraits<util::Vector3D<site_t> >::RegisterMpiDataType()
    {
      return MpiDataTypeTraits<std::array<site_t, 3>>::RegisterMpiDataType();
    }
    template<>
    MPI_Datatype MpiDataTypeTraits<util::Vector3D<U16> >::RegisterMpiDataType()
    {
      return MpiDataTypeTraits<std::array<U16, 3>>::RegisterMpiDataType();
    }

    template<>
    MPI_Datatype MpiDataTypeTraits<util::Vector3D<distribn_t> >::RegisterMpiDataType()
    {
      return MpiDataTypeTraits<std::array<distribn_t, 3>>::RegisterMpiDataType();
    }
}
namespace hemelb::util
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
