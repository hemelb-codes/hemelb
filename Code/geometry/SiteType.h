//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

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
