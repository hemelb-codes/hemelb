//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_FLOWEXTENSION_H
#define HEMELB_REDBLOOD_FLOWEXTENSION_H

#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {

    //! Cylindrical extension at an in/outlet of the vascular system
    class FlowExtension
    {
      public:
        //! \brief Direction of the axis
        //! \details Should be normalized
        util::Vector3D<Dimensionless> normal;
        //! \brief Point on the axis giving the origin
        //! \details Whether this is the left-most or right-most point on the axis depends on the
        //! direction of the normal. In practice, it should be at the intersection between the
        //! cylinder axis and the surface separating the extension from the real flow.
        LatticePosition origin;
        //! Length of the cylinder
        LatticeDistance length;
        //! Radius of the cylinder
        LatticeDistance radius;
        //! Distance within which to fade in/out
        LatticeDistance fadeLength;
    };

    //! Checks whether a cell is inside a flow extension
    bool contains(const FlowExtension &, const LatticePosition &);

    //! Linear weight associated with a point in the cylinder
    Dimensionless linearWeight(FlowExtension const&, LatticePosition const&);

  } // namespace hemelb::redblood
} // namespace hemelb

#endif
