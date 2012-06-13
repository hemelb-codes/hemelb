#ifndef HEMELB_UNITTESTS_GEOMETRY_LATTICEDATATESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_LATTICEDATATESTS_H

#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {

        using namespace hemelb::geometry;
        class NeighbouringLatticeDataTests : public FourCubeBasedTestFixture
        {
            CPPUNIT_TEST_SUITE (NeighbouringLatticeDataTests);
            CPPUNIT_TEST (TestConstruct);
            CPPUNIT_TEST (TestConvertGlobalId);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringLatticeDataTests()
            {
            }

            void setUp()
            {
              FourCubeBasedTestFixture::setUp();

            }

            void tearDown()
            {
              FourCubeBasedTestFixture::tearDown();

            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

            void TestConvertGlobalId()
            {
              // not really a very good test to use a one-proc geometry
              // we need to create a sixteen-cube lattice data test fixture to test this kind of thing properly
              // but we can at least smoke test things with four cube.
              util::Vector3D<site_t> exampleCoord(1,2,3);
              site_t id = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(exampleCoord);
              CPPUNIT_ASSERT_EQUAL(id,static_cast<site_t>(1*16+2*4+3));
              util::Vector3D<site_t> resultCoord;
              latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(id,resultCoord);
              CPPUNIT_ASSERT_EQUAL(resultCoord,exampleCoord);
            }

          private:
        };
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringLatticeDataTests);
      }
    }
  }

#endif
