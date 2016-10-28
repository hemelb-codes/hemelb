
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "unittests/helpers/MockCommsHelper.h"
#include "comm/MpiEnvironment.h"
#include "comm/AsyncConcern.h"

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
                                             public MockCommsHelper
        {
            CPPUNIT_TEST_SUITE ( NeighbouringDataManagerTests);
            CPPUNIT_TEST ( TestConstruct);
            CPPUNIT_TEST ( TestRegisterNeedOneProc);
            CPPUNIT_TEST ( TestShareNeedsOneProc);
            CPPUNIT_TEST ( TestShareConstantDataOneProc);
            CPPUNIT_TEST ( TestShareFieldDataOneProc);
            CPPUNIT_TEST ( TestShareFieldDataOneProcViaIterableAction);

            CPPUNIT_TEST_SUITE_END();

          public:
            NeighbouringDataManagerTests() :
              data(NULL)
            {
            }

            void setUp()
            {
              FourCubeBasedTestFixture::setUp();
              data = &latDat->GetNeighbouringData();
              MockCommsHelper::setUp(1, 0);
              manager = new NeighbouringDataManager(*latDat, *data, commQ);
            }

            void UseRealCommunicator()
            {
              // delete manager;

              // communicator = comm::MpiEnvironment::World();
	      // commQ = comm::Async::New(communicator);
              // manager = new NeighbouringDataManager(*latDat, *data, commQ);
            }

            void tearDown()
            {
              delete manager;
              FourCubeBasedTestFixture::tearDown();
              MockCommsHelper::tearDown();
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

            void TestRegisterNeedOneProc()
            {
              // We imagine that unbeknownst to us, there is a site at x=1,y=1,z=7 which for some unexplained reason we need to know about.
              // The client code has a duty to determine the global index for the site, which it can do,
              // using LatticeData::GetGlobalNoncontiguousSiteIdFromGlobalCoords
              // and LatticeData::GetGlobalCoords
              // for our four cube test case, we only have one proc, of course, so we pretend to require a site from our own proc via this mechanism
              // just to verify it.

              util::Vector3D<site_t> exampleBlockLocalCoord(1, 2, 3);
              util::Vector3D<site_t> exampleBlockCoord(0, 0, 0);
              util::Vector3D<site_t> globalCoord = latDat->GetGlobalCoords(exampleBlockCoord,
                                                                           exampleBlockLocalCoord);
              // to illustrate a typical use, we now add a displacement to the global coords of the local site, which in a real example
              // would take it off-proc
              util::Vector3D<site_t> desiredGlobalCoord = globalCoord + util::Vector3D<site_t>(1,
                                                                                               0,
                                                                                               0);
              // desiredGC = (2,2,3) => id = (2*block + 2) * block + 3
              site_t desiredId =
                  latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(desiredGlobalCoord);
              CPPUNIT_ASSERT_EQUAL(site_t(87), desiredId);

              manager->RegisterNeededSite(desiredId);
              CPPUNIT_ASSERT_EQUAL(static_cast<std::vector<site_t>::size_type> (1),
                                   manager->GetNeededSites().size());
              CPPUNIT_ASSERT_EQUAL(static_cast<site_t> (87), manager->GetNeededSites()[0]);
            }

            void TestShareNeedsOneProc()
            {
	      std::vector<site_t> neededsites = {43};
              manager->RegisterNeededSite(neededsites[0]);
	      
	      // The MapA2A needs a barrier but as the only needs are
	      // local this is true immediately
	      MockComms()->RequireIbarrier([]() { return true; });

	      // We then send this need to ourself
	      MockComms()->RequireSend(neededsites, 0, 0);
	      MockComms()->RequireRecv(neededsites, 0, 0);
              manager->ShareNeeds();
	      MockComms()->ExpectationsAllCompleted();
              CPPUNIT_ASSERT_EQUAL(manager->GetNeedsForProc(0).size(),
                                   static_cast<std::vector<int>::size_type> (1));
              CPPUNIT_ASSERT_EQUAL(manager->GetNeedsForProc(0).front(), static_cast<site_t> (43));
            }

            void TestShareConstantDataOneProc()
            {
	      std::vector<site_t> neededsites = {43};
              manager->RegisterNeededSite(neededsites[0]);
	      
	      MockComms()->RequireIbarrier([]() { return true; });
	      MockComms()->RequireSend(neededsites, 0, 0);
	      MockComms()->RequireRecv(neededsites, 0, 0);
              manager->ShareNeeds();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite
                  = latDat->GetSite(latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(43));
              // It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

              // We should send/receive the site data
              SiteData expectedData = exampleSite.GetSiteData();
              SiteData fixtureData = exampleSite.GetSiteData();

	      MockComms()->RequireSend(expectedData.GetWallIntersectionData(),
				       0, 0);
              MockComms()->RequireRecv(fixtureData.GetWallIntersectionData(),
				       0, 0);
              MockComms()->RequireSend(&expectedData.GetIoletIntersectionData(),
				       1,
				       0,
				       0);
              MockComms()->RequireRecv(&fixtureData.GetIoletIntersectionData(),
				       1,
				       0,
				       0);
              MockComms()->RequireSend(&expectedData.GetIoletId(), 1, 0, 0);
              MockComms()->RequireRecv(&fixtureData.GetIoletId(), 1, 0, 0);
              MockComms()->RequireSend(&expectedData.GetSiteType(), 1, 0, 0);
              MockComms()->RequireRecv(&fixtureData.GetSiteType(), 1, 0, 0);

              MockComms()->RequireSend(exampleSite.GetWallDistances(),
                                   lb::lattices::D3Q15::NUMVECTORS - 1,
                                   0,
                                   0);
              MockComms()->RequireRecv(exampleSite.GetWallDistances(),
                                      lb::lattices::D3Q15::NUMVECTORS - 1,
                                      0,
                                      0);
              MockComms()->RequireSend(&exampleSite.GetWallNormal(), 1, 0, 0);
              MockComms()->RequireRecv(&exampleSite.GetWallNormal(), 1, 0, 0);
              manager->TransferNonFieldDependentInformation();
	      
              MockComms()->ExpectationsAllCompleted();
	      
              NeighbouringSite transferredSite = data->GetSite(43);
              CPPUNIT_ASSERT_EQUAL(exampleSite.GetSiteData(), transferredSite.GetSiteData());
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(exampleSite.GetWallDistances()[direction],
                                     transferredSite.GetWallDistances()[direction]);
              }
              CPPUNIT_ASSERT_EQUAL(exampleSite.GetWallNormal(), transferredSite.GetWallNormal());
            }

            void TestShareFieldDataOneProc()
            {
              

              site_t targetGlobalOneDIdx = 43;
              LatticeVector targetGlobalThreeDIdx = latDat->GetSiteCoordsFromSiteId(targetGlobalOneDIdx);
              site_t targetLocalIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(targetGlobalOneDIdx);

              for (unsigned int direction = 0; direction < 3; direction++)
                CPPUNIT_ASSERT_EQUAL(site_t(1), targetGlobalThreeDIdx[direction]);
	      
	      MockComms()->RequireIbarrier([]() { return true; });
	      std::vector<site_t> neededsites = {targetGlobalOneDIdx};
	      MockComms()->RequireSend(neededsites, 0, 0);
	      MockComms()->RequireRecv(neededsites, 0, 0);

              manager->RegisterNeededSite(targetGlobalOneDIdx);
              manager->ShareNeeds();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);

              // It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

              MockComms()->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
				       lb::lattices::D3Q15::NUMVECTORS,
				       0,
				       0);

              std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
              MockComms()->RequireRecv(&(receivedFOld[0]),
				       lb::lattices::D3Q15::NUMVECTORS,
				       0,
				       0);

              manager->TransferFieldDependentInformation();
	      MockComms()->ExpectationsAllCompleted();

              NeighbouringSite transferredSite = data->GetSite(targetGlobalOneDIdx);
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(receivedFOld[direction],
                                     transferredSite.GetFOld<lb::lattices::D3Q15> ()[direction]);
              }
            }

            void TestShareFieldDataOneProcViaIterableAction()
            {
              site_t targetGlobalOneDIdx = 43;
              site_t targetLocalIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(targetGlobalOneDIdx);
	      
	      MockComms()->RequireIbarrier([]() { return true; });
	      std::vector<site_t> neededsites = {targetGlobalOneDIdx};
	      MockComms()->RequireSend(neededsites, 0, 0);
	      MockComms()->RequireRecv(neededsites, 0, 0);

              manager->RegisterNeededSite(targetGlobalOneDIdx);
              manager->ShareNeeds();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);
              // It should arrive in the NeighbouringDataManager, from
              // the values sent from the localLatticeData

              MockComms()->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
				       lb::lattices::D3Q15::NUMVECTORS,
				       0,
				       0);
              std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
              MockComms()->RequireRecv(&(receivedFOld[0]),
				       lb::lattices::D3Q15::NUMVECTORS,
				       0,
				       0);
	      
              net::phased::StepManager stepManager;
              stepManager.RegisterIteratedActorSteps(*manager, 0);
	      comm::AsyncConcern netConcern = comm::AsyncConcern(commQ);
              stepManager.RegisterCommsForAllPhases(netConcern);
              stepManager.CallActions();
              MockComms()->ExpectationsAllCompleted();

              NeighbouringSite transferredSite = data->GetSite(targetGlobalOneDIdx);
              for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
              {
                CPPUNIT_ASSERT_EQUAL(receivedFOld[direction],
                                     transferredSite.GetFOld<lb::lattices::D3Q15> ()[direction]);
              }
            }

          private:
            NeighbouringDataManager *manager;
            NeighbouringLatticeData *data;

        };
        // CPPUNIT_EXTRA_LINE
        CPPUNIT_TEST_SUITE_REGISTRATION ( NeighbouringDataManagerTests);
      }
    }
  }
}

#endif
