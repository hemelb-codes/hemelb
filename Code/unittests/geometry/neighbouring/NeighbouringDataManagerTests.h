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
              MockNetHelper::setUp(1,0);
              data = new NeighbouringLatticeData(latDat->GetLatticeInfo());
              manager = new NeighbouringDataManager(*latDat, *data, *netMock);
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
              // The client code has a duty to determine the global index for the site, which it can do,
              // using LatticeData::GetGlobalNoncontiguousSiteIdFromGlobalCoords
              // and LatticeData::GetGlobalCoords
              // for our four cube test case, we only have one proc, of course, so we pretend to require a site from our own proc via this mechanism
              // just to verify it.

              util::Vector3D<site_t> exampleBlockLocalCoord(1, 2, 3);
              util::Vector3D<site_t> exampleBlockCoord(0, 0, 0);
              util::Vector3D<site_t> globalCoord=latDat->GetGlobalCoords(exampleBlockCoord,exampleBlockLocalCoord);
              // to illustrate a typical use, we now add a displacement to the global coords of the local site, which in a real example
              // would take it off-proc
              util::Vector3D<site_t> desiredGlobalCoord=globalCoord+util::Vector3D<site_t>(1,0,0);
              site_t desiredId = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(desiredGlobalCoord);
              CPPUNIT_ASSERT_EQUAL(desiredId,static_cast<site_t>(43)); // 43 = 2*16+2*4+3

              manager->RegisterNeededSite(desiredId);
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
