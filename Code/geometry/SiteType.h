// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITETYPE_H
#define HEMELB_GEOMETRY_SITETYPE_H

namespace hemelb::geometry
{
    /*
     * Note: The implementation contains the MPI type registration code.
     */
    enum SiteType: unsigned
    {
      // These must be consistent with the setup tool
      SOLID_TYPE = 0U,
      FLUID_TYPE = 1U,
      INLET_TYPE = 2U,
      OUTLET_TYPE = 3U
    };

}

#ifdef HEMELB_CODE

// If we are building the main app, let MPI use the underlying type of
// the enum for comms.
#include "net/MpiDataType.h"

namespace hemelb::geometry {
  inline MPI_Datatype MpiDataType(SiteType const&) {
    return net::MpiDataType<unsigned>();
  }
}

#endif

#endif // HEMELB_GEOMETRY_SITETYPE_H
