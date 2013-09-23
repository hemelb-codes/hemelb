//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/mpi.h"
#include "geometry/SiteType.h"

namespace hemelb
{
  namespace net {
  template<>
  MPI_Datatype MpiDataTypeTraits<geometry::SiteType>::RegisterMpiDataType()
  {
    return MpiDataTypeTraits<int>::RegisterMpiDataType();
  }
  }
}
