#ifndef HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    /**
     * Class for selecting the whole geometry.
     */
    class WholeGeometrySelector : public GeometrySelector
    {
      protected:
        /**
         * Returns true for all locations.
         *
         * @param data
         * @param location
         * @return
         */
        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);
    };
  }
}

#endif /* HEMELB_EXTRACTION_WHOLEGEOMETRYSELECTOR_H */
