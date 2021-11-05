// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GEOMETRYSITELINK_H
#define HEMELB_GEOMETRY_GEOMETRYSITELINK_H

#include "io/formats/geometry.h"

namespace hemelb
{
  namespace geometry
  {
    //! Model of the data read in about a link from one site in one direction.
    struct GeometrySiteLink
    {
      public:
        io::formats::geometry::CutType type;

        //! Default constructor. Has no intersection, nonsense values for intersection distance
        //! and iolet id.
        GeometrySiteLink() :
            type(io::formats::geometry::CutType::NONE), distanceToIntersection(-1.0), ioletId(-1)
        {
        }

        //! The proportion of the lattice vector traversed before an intersection is hit.
        float distanceToIntersection;

        //! The identifier of the inlet or outlet hit along the lattice vector (if one is hit).
        int ioletId;
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYSITELINK_H */
