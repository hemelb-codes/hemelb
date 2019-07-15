// This file is part of HemeLB and is Copyright (C)
// // the HemeLB team and/or their institutions, as detailed in the
// // file AUTHORS. This software is provided under the terms of the
// // license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"
#include "geometry/SiteType.h"

namespace hemelb
{
  namespace extraction
  {
    class IoletGeometrySelector : public GeometrySelector
    {
      public:

        IoletGeometrySelector(int ioletId, geometry::SiteType ioletType);

        int GetIoletId() const;

        geometry::SiteType GetIoletType() const;

      protected:

        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:

        const int ioletId;
        const geometry::SiteType ioletType;
    };
  }
}

#endif /* HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H */
