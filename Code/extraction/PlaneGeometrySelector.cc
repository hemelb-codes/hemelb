
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PlaneGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    PlaneGeometrySelector::PlaneGeometrySelector(const util::Vector3D<float>& point,
                                                 const util::Vector3D<float>& normal) :
      planePoint(point), normal(normal.GetNormalised()), radius(0.)
    {

    }

    /**
     * Constructor makes a plane geometry object with given normal, about a given
     * point, with given radius.
     * @param point
     * @param normal
     * @param radius
     */
    PlaneGeometrySelector::PlaneGeometrySelector(const util::Vector3D<float>& point,
                                                 const util::Vector3D<float>& normal,
                                                 float radius) :
      planePoint(point), normal(normal.GetNormalised()), radius(radius)
    {

    }

    const util::Vector3D<float>& PlaneGeometrySelector::GetPoint() const
    {
      return planePoint;
    }

    const util::Vector3D<float>& PlaneGeometrySelector::GetNormal() const
    {
      return normal;
    }

    float PlaneGeometrySelector::GetRadius() const
    {
      return radius;
    }

    bool PlaneGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                 const util::Vector3D<site_t>& location)
    {
      util::Vector3D<float> coords = util::Vector3D<float>(location) * data.GetVoxelSize() + data.GetOrigin();

      const float perpendicularDistance = (coords - planePoint).Dot(normal);

      if (std::abs(perpendicularDistance) > (0.5 * data.GetVoxelSize()))
      {
        return false;
      }

      // Return true if using infinite radius of the plane.
      if (radius <= 0.)
      {
        return true;
      }

      const float radiusOfPointSquared =
          ( (coords - normal * perpendicularDistance) - planePoint).GetMagnitudeSquared();

      return radiusOfPointSquared <= radius * radius;
    }
  }
}
