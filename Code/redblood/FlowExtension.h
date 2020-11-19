// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_FLOWEXTENSION_H
#define HEMELB_REDBLOOD_FLOWEXTENSION_H

#include <memory>
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {
    struct Cylinder
    {
        //! \brief Direction of the axis
        //! \details Should be normalized
        LatticePosition normal;
        //! \brief Point on the axis giving the origin
        //! \details Whether this is the left-most or right-most point on the axis depends on the
        //! direction of the normal. In practice, it should be at the intersection between the
        //! cylinder axis and the surface separating the extension from the real flow.
        LatticePosition origin;
        //! Radius of the cylinder
        LatticeDistance radius;
        //! Length of the cylinder
        LatticeDistance length;
    };

    //! Cylindrical extension at an in/outlet of the vascular system
    class FlowExtension : public Cylinder
    {
      public:
        //! Distance within which to fade in/out
        LatticeDistance fadeLength;
        FlowExtension(LatticePosition const &n0, LatticePosition const &gamma, LatticeDistance l,
                      LatticeDistance r, LatticeDistance fl) :
            Cylinder( { n0, gamma, r, l }), fadeLength(fl)
        {
        }
        FlowExtension() :
            Cylinder( { LatticePosition(1, 0, 0), LatticePosition(0, 0, 0), 1, 0 }), fadeLength(1)
        {
        }
    };

    //! Checks whether a cell is inside a flow extension
    bool contains(const Cylinder &, const LatticePosition &);
    //! Checks whether a cell is inside a flow extension
    inline bool contains(std::shared_ptr<Cylinder const> cylinder, const LatticePosition &position)
    {
      return contains(*cylinder, position);
    }

    //! Linear weight associated with a point in the cylinder
    Dimensionless linearWeight(FlowExtension const&, LatticePosition const&);

  } // namespace hemelb::redblood
} // namespace hemelb

#endif
