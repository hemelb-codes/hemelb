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
