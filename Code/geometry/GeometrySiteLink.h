
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GEOMETRYSITELINK_H
#define HEMELB_GEOMETRY_GEOMETRYSITELINK_H

namespace hemelb
{
  namespace geometry
  {
    //! Model of the data read in about a link from one site in one direction.
    struct GeometrySiteLink
    {
      public:
        //! Enumeration of the different intersections the link might make between the current
        //! site and the next lattice point in this direction: no intersection,
        //! intersection with a vessel wall and intersection with an inlet or outlet.
        //! @todo #598 remove enumerated type and replace with equivalent in io::formats::geometry
        enum IntersectionType
        {
          NO_INTERSECTION = 0,
          WALL_INTERSECTION = 1,
          INLET_INTERSECTION = 2,
          OUTLET_INTERSECTION = 3
        } type;

        //! Default constructor. Has no intersection, nonsense values for intersection distance
        //! and iolet id.
        GeometrySiteLink() :
            type(NO_INTERSECTION), distanceToIntersection(-1.0), ioletId(-1)
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
