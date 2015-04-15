// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
            && (x.GetIoletIntersectionData() == y.GetIoletIntersectionData())
            && (x.GetSiteType() == y.GetSiteType()) && (x.GetIoletId() == y.GetIoletId());
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
