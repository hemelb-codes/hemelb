// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_EQUALITYSITEDATA_H
#define HEMELB_TESTS_HELPERS_EQUALITYSITEDATA_H

#include "geometry/SiteData.h"

namespace hemelb
{
  namespace geometry {
    inline bool operator==(const hemelb::geometry::SiteData& x, const hemelb::geometry::SiteData& y) {
      return (x.GetWallIntersectionData() == y.GetWallIntersectionData()) && (x.GetIoletIntersectionData() == y.GetIoletIntersectionData()) && (x.GetSiteType()
            == y.GetSiteType()) && (x.GetIoletId() == y.GetIoletId());
    }
  }
}

#endif // ONCE
