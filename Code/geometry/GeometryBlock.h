// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
