
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PLANEGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_PLANEGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Selects a geometry that forms a squat cylinder with height 0.5 lattice units
     * around a planar circle with specified normal and centre, and optionally specified
     * radius (assumed to be infinite when absent).
     */
    class PlaneGeometrySelector : public GeometrySelector
    {
      public:
        /**
         * Constructor makes an infinite plane geometry object with given normal, about a given
         * point.
         * @param point
         * @param normal
         */
        PlaneGeometrySelector(const util::Vector3D<float>& point, const util::Vector3D<float>& normal);

        /**
         * Constructor makes a plane geometry object with given normal, about a given
         * point, with given radius.
         * @param point
         * @param normal
         * @param radius
         */
        PlaneGeometrySelector(const util::Vector3D<float>& point, const util::Vector3D<float>& normal, float radius);

        /**
         * Returns a point that lies on the plane.
         * @return
         */
        const util::Vector3D<float>& GetPoint() const;

        /**
         * Returns the plane normal.
         */
        const util::Vector3D<float>& GetNormal() const;

        /**
         * Returns the radius of the plane.
         */
        float GetRadius() const;

      protected:
        /**
         * Returns true for any location within 0.5 lattice units of the plane / squat cylinder.
         *
         * @param data
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:
        /**
         * A point on the plane.
         */
        const util::Vector3D<float> planePoint;

        /**
         * The plane normal.
         */
        const util::Vector3D<float> normal;

        /**
         * The radius around the planePoint to select. Radius <= 0 is taken as infinite.
         */
        const float radius;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PLANEGEOMETRYSELECTOR_H */
