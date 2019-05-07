// This file is part of HemeLB and is Copyright (C)
// // the HemeLB team and/or their institutions, as detailed in the
// // file AUTHORS. This software is provided under the terms of the
// // license in the file LICENSE.

#include "extraction/IoletGeometrySelector.h"

namespace hemelb
{
  namespace extraction
  {
    IoletGeometrySelector::IoletGeometrySelector(int ioletId, unsigned int ioletType) :
      ioletId(ioletId), ioletType(ioletType)
    {

    }

    int IoletGeometrySelector::GetIoletId() const
    {
      return ioletId;
    }

    unsigned int IoletGeometrySelector::GetIoletType() const
    {
      return ioletType;
    }

    bool IoletGeometrySelector::IsWithinGeometry(const extraction::IterableDataSource& data,
                                                 const util::Vector3D<site_t>& location)
    {
       bool result = false;
       result = data.IsIoletSite(location, ioletId, ioletType);

       return result;
    }
  }
}
