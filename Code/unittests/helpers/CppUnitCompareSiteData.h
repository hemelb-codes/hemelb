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
        return (x.GetIntersectionData() == y.GetIntersectionData()) && (x.GetOtherRawData() == y.GetOtherRawData());
      }

      // Note this vector print doesn't visually distinguish between ("1" "2" "3") and("1, 2" "3").
      static std::string toString(const hemelb::geometry::SiteData& value)
      {
        std::stringstream output;
        output << value.GetIntersectionData() << " , " << value.GetOtherRawData() << std::flush;
        return output.str();
      }
  };
}
#endif
