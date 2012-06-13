#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H

#include "geometry/neighbouring/NeighbouringDataManager.h"

namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {
      namespace neighbouring
      {
        using namespace hemelb::geometry::neighbouring;
        class NeighbouringDataManagerTests : public FourCubeBasedTestFixture
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
            }

            void tearDown()
            {
              FourCubeBasedTestFixture::tearDown();
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

          private:
            net::NetMock *netMock;
            topology::Communicator *communicatorMock;
            NeighbouringDataManager *manager;

        };
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringDataManagerTests);
      }
    }
  }
}

#endif
