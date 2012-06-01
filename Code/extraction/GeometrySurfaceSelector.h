#ifndef HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H
#define HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {

    class GeometrySurfaceSelector : public hemelb::extraction::GeometrySelector
    {
      protected:

        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);
    };

  } /* namespace extraction */
} /* namespace hemelb */
#endif /* HEMELB_EXTRACTION_GEOMETRYSURFACESELECTOR_H */
