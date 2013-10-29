//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_GEOMETRY_GEOMETRYPEEKER_H
#define HEMELB_GEOMETRY_GEOMETRYPEEKER_H

#include "geometry/GeometryReaderBase.h"

namespace hemelb
{
  namespace geometry
  {
    class GeometryPeeker: public GeometryReaderBase
    {
      public:
        GeometryPeeker(const std::string& gmyFileName);
        inline const LatticeDistance& GetVoxelSize()
        {
          return voxelSize;
        }
        inline const LatticePosition& GetOrigin()
        {
          return origin;
        }

      private:
        LatticeDistance voxelSize;
        LatticePosition origin;
    };
  }
}
#endif // HEMELB_GEOMETRY_GEOMETRYPEEKER_H
