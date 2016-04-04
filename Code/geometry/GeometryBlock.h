
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GEOMETRYBLOCK_H
#define HEMELB_GEOMETRY_GEOMETRYBLOCK_H

#include <vector>
#include "geometry/GeometrySite.h"

namespace hemelb
{
  namespace geometry
  {
    /***
     * Model of the information stored for a block in a geometry file.
     * Just gives the array of sites
     */
    struct BlockReadResult
    {
      public:
        std::vector<GeometrySite> Sites;
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYBLOCK_H */
