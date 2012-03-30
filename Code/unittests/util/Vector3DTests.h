#ifndef HEMELB_UNITTESTS_UTIL_VECTOR3DTESTS_H
#define HEMELB_UNITTESTS_UTIL_VECTOR3DTESTS_H

#include <cppunit/TestFixture.h>
#include "util/Vector3D.h"

#include <boost/type_traits.hpp>

namespace hemelb
{
  namespace unittests
  {
    namespace util
    {
      class Vector3DTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(Vector3DTests);
          CPPUNIT_TEST(TestCastsInVector3DProduct);CPPUNIT_TEST_SUITE_END();
        public:
          void TestCastsInVector3DProduct()
          {
            const double dblMax = std::numeric_limits<double>::max();
            const unsigned uintMax = std::numeric_limits<unsigned>::max();

            {
              // (int, float) -> float
              hemelb::util::Vector3D<int> foo(-1, 0, 1);
              float bar = 0.8;
              hemelb::util::Vector3D<float> baz = foo * bar;

              CPPUNIT_ASSERT_EQUAL(-0.8f, baz[0]);
              CPPUNIT_ASSERT_EQUAL(0.0f, baz[1]);
              CPPUNIT_ASSERT_EQUAL(0.8f, baz[2]);
            }
            {
              // (float, double) -> double
              hemelb::util::Vector3D<float> foo(0.0f, 1.0f, -1.0f);
              double bar = dblMax;
              hemelb::util::Vector3D<double> baz = foo * bar;

              CPPUNIT_ASSERT_EQUAL(0.0, baz[0]);
              CPPUNIT_ASSERT_EQUAL(dblMax, baz[1]);
              CPPUNIT_ASSERT_EQUAL(-dblMax, baz[2]);
            }
            {
              // (int, unsigned) -> unsigned
              hemelb::util::Vector3D<int> foo(0, 2, 2);
              unsigned bar = uintMax / 2u;
              hemelb::util::Vector3D<unsigned> baz = foo * bar;

              CPPUNIT_ASSERT_EQUAL(0u, baz[0]);
              CPPUNIT_ASSERT_EQUAL(uintMax, baz[1] + uintMax % 2);
              CPPUNIT_ASSERT_EQUAL(uintMax, baz[2] + uintMax % 2);
            }
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(Vector3DTests);

    }
  }
}
#endif /* HEMELB_UNITTESTS_UTIL_VECTOR3DTESTS_H */
