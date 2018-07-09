
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_EXTRACTION_SURFACEPOINTSELECTOR_H
#define HEMELB_EXTRACTION_SURFACEPOINTSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {

    /**
     * This class defines a geometry selector for a point in the geometry surface. By geometry surface we mean the
     * original surface passed into the setup tool (a .stl file possibly). Since such point is unlikely to map directly
     * onto a lattice site, the class will select any lattice site known to be a wall (in LatticeData terms) and within
     * a sqrt(3) times voxel size radius. This can obviously lead to multiple lattice sites being selected. The user will
     * have to handle this at postprocessing time.
     */
    class SurfacePointSelector : public GeometrySelector
    {
      public:
        /**
         * Constructor takes a point in the surface of the original .stl file
         * @param surfacePoint
         */
        SurfacePointSelector(const util::Vector3D<float>& surfacePoint);

        /**
         * Get the surface point.
         * @return surface point
         */
        const util::Vector3D<float>& GetPoint() const;

      protected:
        /**
         * Returns true if the given location is within the selection.
         * @param data data source abstraction
         * @param location lattice site coordinates to evaluate the selector on
         * @return whether location is within the selection
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:
        /** Coordinates of the surface point to be selected. */
        const util::Vector3D<float> surfacePoint;
        /** Maximum distance between two sites in lattice units. */
        static const float maxDistanceBetweenSitesLatticeUnits;
    };

  } /* namespace extraction */
} /* namespace hemelb */
#endif /* HEMELB_EXTRACTION_SURFACEPOINTSELECTOR_H */
