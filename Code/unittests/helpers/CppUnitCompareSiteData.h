
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_CPPUNITCOMPARESITEDATA_H
#define HEMELB_UNITTESTS_HELPERS_CPPUNITCOMPARESITEDATA_H
#include <cppunit/TestFixture.h>
#include "geometry/SiteData.h"

namespace CPPUNIT_NS
{
  template<> struct assertion_traits<hemelb::geometry::SiteData>
  {
      static bool equal(const hemelb::geometry::SiteData& x, const hemelb::geometry::SiteData& y)
      {
        return (x.GetWallIntersectionData() == y.GetWallIntersectionData())
            && (x.GetIoletIntersectionData() == y.GetIoletIntersectionData()) && (x.GetSiteType()
            == y.GetSiteType()) && (x.GetIoletId() == y.GetIoletId());
      }

      // Note this vector print doesn't visually distinguish between ("1" "2" "3") and("1, 2" "3").
      static std::string toString(const hemelb::geometry::SiteData& value)
      {
        std::stringstream output;
        output << value.GetSiteType() << ", " << value.GetWallIntersectionData() << ", "
            << value.GetIoletIntersectionData() << ", " << value.GetIoletId() << std::flush;
        return output.str();
      }
  };
}
#endif
