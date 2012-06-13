#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H

#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "unittests/helpers/MockNetHelper.h"

namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class NeighbouringDataManagerTests : public FourCubeBasedTestFixture, public MockNetHelper
        {
            CPPUNIT_TEST_SUITE (NeighbouringDataManagerTests);
            CPPUNIT_TEST (TestConstruct);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringDataManagerTests()
            {
            }

            void setUp()
            {
              FourCubeBasedTestFixture::setUp();
              manager=new NeighbouringDataManager();
            }

            void tearDown()
            {
              FourCubeBasedTestFixture::tearDown();
              MockNetHelper::tearDown();
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

          private:
            NeighbouringDataManager *manager;

        };
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringDataManagerTests);
      }
    }
  }
}

#endif
