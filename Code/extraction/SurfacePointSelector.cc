#include "SurfacePointSelector.h"

#include <iostream>

namespace hemelb
{
  namespace extraction
  {
    SurfacePointSelector::SurfacePointSelector(const util::Vector3D<float>& surfacePoint) :
        surfacePoint(surfacePoint), maxDistanceBetweenSitesLatticeUnits(sqrt(3.0))
    {
    }

    util::Vector3D<float> SurfacePointSelector::GetPoint() const
    {
      return surfacePoint;
    }

    bool SurfacePointSelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                const util::Vector3D<site_t>& location)
    {
      // Don't bother checking distance if current site is not marked as edge
      if (!data.IsEdgeSite(location))
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
