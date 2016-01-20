
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/StraightLineGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    StraightLineGeometrySelector::StraightLineGeometrySelector(const util::Vector3D<float>& endpoint1,
                                               const util::Vector3D<float>& endpoint2) :
      endpoint1(endpoint1), lineVector(endpoint2 - endpoint1), lineLength(lineVector.GetMagnitude())
    {

    }

    util::Vector3D<float> StraightLineGeometrySelector::GetEndpoint1() const
    {
      return endpoint1;
    }

    util::Vector3D<float> StraightLineGeometrySelector::GetEndpoint2() const
    {
      return lineVector + endpoint1;
    }

    bool StraightLineGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                        const util::Vector3D<site_t>& location)
    {
      util::Vector3D<float> coords = util::Vector3D<float>(location) * data.GetVoxelSize() + data.GetOrigin();

      const float lengthAlongLine = (lineVector.Dot(coords - endpoint1)) / lineLength;

      if (lengthAlongLine < 0. || lengthAlongLine > lineLength)
      {
        return false;
      }

      // Use magnitude squared as it saves a sqrt operation.
      const float perpendicularDistanceSquared =
          ( (endpoint1 + lineVector * lengthAlongLine / lineLength) - coords).GetMagnitudeSquared();

      // This is chosen so that a line as far as possible from all lattice points will
      // still gather some data.
      // The worst case is when the line is lattice-aligned and is passing through the centre
      // of voxels. The perpendicular distance from the line here is the distance from the corner
      // of a voxel-sized square to its centre.
      return perpendicularDistanceSquared <= (2.0 * 0.5 * 0.5 * data.GetVoxelSize() * data.GetVoxelSize());
    }
  }
}
