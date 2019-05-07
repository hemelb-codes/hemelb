// This file is part of HemeLB and is Copyright (C)
// // the HemeLB team and/or their institutions, as detailed in the
// // file AUTHORS. This software is provided under the terms of the
// // license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H
#define HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H

#include "extraction/GeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    class IoletGeometrySelector : public GeometrySelector
    {
      public:

        IoletGeometrySelector(int ioletId, unsigned int ioletType);

        int GetIoletId() const;

        unsigned int GetIoletType() const;

      protected:

        bool IsWithinGeometry(const extraction::IterableDataSource& data, const util::Vector3D<site_t>& location);

      private:

        const int ioletId;
        const unsigned ioletType;
    };
  }
}

#endif /* HEMELB_EXTRACTION_IOLETGEOMETRYSELECTOR_H */
