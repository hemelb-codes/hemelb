//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/FlowExtension.h"
#include "redblood/Mesh.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {

    //! Checks whether a cell is inside a flow extension
    bool contains(const Cylinder & flowExt, const LatticePosition& point)
    {
      assert(std::abs(flowExt.normal.GetMagnitude() - 1e0) < 1e-8);
      assert(flowExt.length > 1e-12);
      assert(flowExt.radius > 1e-12);
      // Vector from centre of start of cylinder to point
      LatticePosition const pd = point - flowExt.origin;

      // distance from pos to point
      LatticeDistance const dot = pd.Dot(flowExt.normal);

      // If the (squared) distance is less than 0 then the point is behind the
      // cylinder cap at pos.
      // If the (squared) distance is greater than the squared length of the
      // cylinder then the point is beyond the end of the cylinder.
      if (dot < 0.0 || dot > flowExt.length)
        return false;

      // Point lies between the parallel caps so calculate the (squared) distance
      // between the point and the line between a and b
      LatticeDistance const dist = pd.GetMagnitudeSquared() - dot * dot;
      return (dist <= flowExt.radius * flowExt.radius);
    }

    Dimensionless linearWeight(FlowExtension const& flowExt, LatticePosition const& position)
    {
      assert(std::abs(flowExt.normal.GetMagnitude() - 1e0) < 1e-8);
      assert(flowExt.length > 1e-12);
      assert(flowExt.radius > 1e-12);

      if (not contains(flowExt, position))
        return 0e0;

      auto const z = flowExt.normal.Dot(position - flowExt.origin);
      return std::max(0e0, 1e0 - z / flowExt.fadeLength);
    }
  } // namespace hemelb::redblood
} // namespace hemelb
