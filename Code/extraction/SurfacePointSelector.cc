
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "extraction/SurfacePointSelector.h"

namespace hemelb
{
  namespace extraction
  {
    const float SurfacePointSelector::maxDistanceBetweenSitesLatticeUnits = sqrt(3.0);

    SurfacePointSelector::SurfacePointSelector(const util::Vector3D<float>& surfacePoint) :
        surfacePoint(surfacePoint)
    {
    }

    const util::Vector3D<float>& SurfacePointSelector::GetPoint() const
    {
      return surfacePoint;
    }

    bool SurfacePointSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                const util::Vector3D<site_t>& location)
    {
      // Don't bother checking distance if current site is not marked as wall
      if (!data.IsWallSite(location))
      {
        return false;
      }

      util::Vector3D<float> coords = util::Vector3D<float>(location) * data.GetVoxelSize() + data.GetOrigin();
      float distanceToSurfacePointLatticeUnits = (coords - surfacePoint).GetMagnitude() / data.GetVoxelSize();

      // Ideally, we would like to check whether the distance is strictly smaller than sqrt(3.0). Unfortunately
      // this is not possible since there's too much floating point arithmetic to yield an exact result and the
      // behaviour could become undefined.
      return (distanceToSurfacePointLatticeUnits <= maxDistanceBetweenSitesLatticeUnits);
    }

  } /* namespace extraction */
} /* namespace hemelb */
