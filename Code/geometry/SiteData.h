
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITEDATA_H
#define HEMELB_GEOMETRY_SITEDATA_H

#include "geometry/SiteDataBare.h"
#include "net/MpiDataType.h"

namespace hemelb
{
  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<geometry::SiteType>::RegisterMpiDataType();
  }
}

#endif /* HEMELB_GEOMETRY_SITEDATA_H */
