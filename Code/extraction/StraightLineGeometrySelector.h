
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_STRAIGHTLINEGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_STRAIGHTLINEGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Selects a geometry that forms a thin cylinder with radius 0.5 lattice units
     * around a line segment between two given points.
     */
    class StraightLineGeometrySelector : public GeometrySelector
    {
      public:
        /**
         * Constructor makes a line geometry object with two endpoints.
         * @param endpoint1
         * @param endpoint2
         */
        StraightLineGeometrySelector(const util::Vector3D<float>& endpoint1, const util::Vector3D<float>& endpoint2);

        /**
         * Get the first endpoint of the line.
         * @return
         */
        util::Vector3D<float> GetEndpoint1() const;

        /**
         * Get the other endpoint of the line.
         * @return
         */
        util::Vector3D<float> GetEndpoint2() const;

      protected:
        /**
         * Returns true for any location within 0.5 lattice units of the line.
         *
         * @param latticeData
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:
        /**
         * The initial point of the line and the line vector.
         */
        const util::Vector3D<float> endpoint1, lineVector;
        /**
         * The length of the line from endpoint1 to endpoint2.
         */
        const float lineLength;
    };
  }
}

#endif /* HEMELB_EXTRACTION_STRAIGHTLINEGEOMETRYSELECTOR_H */
