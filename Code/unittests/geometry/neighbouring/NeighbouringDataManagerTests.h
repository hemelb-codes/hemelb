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
        class NeighbouringDataManagerTests : public FourCubeBasedTestFixture,
                                             public MockNetHelper
        {
            CPPUNIT_TEST_SUITE (NeighbouringDataManagerTests);
            CPPUNIT_TEST (TestConstruct);
            CPPUNIT_TEST (TestRegisterNeed);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringDataManagerTests() :
                data(NULL)
            {
            }

            void setUp()
            {
              FourCubeBasedTestFixture::setUp();
              data = new NeighbouringLatticeData(latDat->GetLatticeInfo());
              manager = new NeighbouringDataManager(*latDat, *data);
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

            void TestRegisterNeed()
            {
              // We imagine that unbeknownst to us, there is a site at x=1,y=1,z=7 which for some unexplained reason we need to know about.
              // The client code has a duty to determine the global index for the site, which it can do, in this way:
            }

          private:
            NeighbouringDataManager *manager;
            NeighbouringLatticeData *data;

        };
        CPPUNIT_TEST_SUITE_REGISTRATION (NeighbouringDataManagerTests);
      }
    }
  }
}

#endif
