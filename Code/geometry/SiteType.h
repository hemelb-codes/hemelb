
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITETYPE_H
#define HEMELB_GEOMETRY_SITETYPE_H

namespace hemelb
{
  namespace geometry
  {
    /*
     * Note: The implementation contains the MPI type registration code.
     */
    enum SiteType
    {
      // These must be consistent with the setup tool
      SOLID_TYPE = 0U,
      FLUID_TYPE = 1U,
      INLET_TYPE = 2U,
      OUTLET_TYPE = 3U
    };

  }

}
#endif // HEMELB_GEOMETRY_SITETYPE_H
