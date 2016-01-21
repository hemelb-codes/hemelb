
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#define HEMELB_UNITTESTS_GEOMETRY_NEIGHBOURING_NEIGHBOURINGDATAMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "net/phased/NetConcern.h"
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
              MockNetHelper::setUp(1, 0);
              manager = new NeighbouringDataManager(*latDat, *data, *netMock);
            }

            void tearDown()
            {
              delete manager;
              FourCubeBasedTestFixture::tearDown();
              MockNetHelper::tearDown();
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
              manager->RegisterNeededSite(43);

              // We should receive a signal that we need one from ourselves
              std::vector<int> countOfNeedsToZeroFromZero;
              countOfNeedsToZeroFromZero.push_back(1); // expectation
              std::vector<int> countOfNeedsFromZeroToZero;
              countOfNeedsFromZeroToZero.push_back(1); //fixture
              netMock->RequireSend(&countOfNeedsToZeroFromZero.front(), 1, 0, "CountToSelf");
              netMock->RequireReceive(&countOfNeedsFromZeroToZero.front(), 1, 0, "CountFromSelf");

              // Once we've received the signal that we need one from ourselves, we should receive that one.
              std::vector<site_t> needsShouldBeSentToSelf;
              std::vector<site_t> needsShouldBeReceivedFromSelf;
              needsShouldBeSentToSelf.push_back(43); //expectation
              needsShouldBeReceivedFromSelf.push_back(43); //fixture
              netMock->RequireSend(&needsShouldBeSentToSelf.front(), 1, 0, "NeedToSelf");
              netMock->RequireReceive(&needsShouldBeSentToSelf.front(), 1, 0, "NeedFromSelf");

              manager->ShareNeeds();
              netMock->ExpectationsAllCompleted();

              CPPUNIT_ASSERT_EQUAL(manager->GetNeedsForProc(0).size(),
                                   static_cast<std::vector<int>::size_type> (1));
              CPPUNIT_ASSERT_EQUAL(manager->GetNeedsForProc(0).front(), static_cast<site_t> (43));
            }

            void TestShareConstantDataOneProc()
            {
              // As for ShareNeeds test, set up the site as needed from itself.
              std::vector<int> countOfNeedsToZeroFromZero;
              countOfNeedsToZeroFromZero.push_back(1); // expectation
              std::vector<int> countOfNeedsFromZeroToZero;
              countOfNeedsFromZeroToZero.push_back(1); //fixture
              netMock->RequireSend(&countOfNeedsToZeroFromZero.front(), 1, 0, "CountToSelf");
              netMock->RequireReceive(&countOfNeedsFromZeroToZero.front(), 1, 0, "CountFromSelf");

              std::vector<site_t> needsShouldBeSentToSelf;
              std::vector<site_t> needsShouldBeReceivedFromSelf;
              needsShouldBeSentToSelf.push_back(43); //expectation
              needsShouldBeReceivedFromSelf.push_back(43); //fixture
              netMock->RequireSend(&needsShouldBeSentToSelf.front(), 1, 0, "NeedToSelf");
              netMock->RequireReceive(&needsShouldBeSentToSelf.front(), 1, 0, "NeedFromSelf");

              manager->RegisterNeededSite(43);
              manager->ShareNeeds();
              netMock->ExpectationsAllCompleted();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite
                  = latDat->GetSite(latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(43));
              // It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

              // We should send/receive the site data
              SiteData expectedData = exampleSite.GetSiteData();
              SiteData fixtureData = exampleSite.GetSiteData();

              netMock->RequireSend(&expectedData.GetWallIntersectionData(),
                                   1,
                                   0,
                                   "WallIntersectionDataToSelf");
              netMock->RequireReceive(&fixtureData.GetWallIntersectionData(),
                                      1,
                                      0,
                                      "WallIntersectionDataFromSelf");
              netMock->RequireSend(&expectedData.GetIoletIntersectionData(),
                                   1,
                                   0,
                                   "IoletIntersectionDataToSelf");
              netMock->RequireReceive(&fixtureData.GetIoletIntersectionData(),
                                      1,
                                      0,
                                      "IoletIntersectionDataFromSelf");
              netMock->RequireSend(&expectedData.GetIoletId(), 1, 0, "IoletIdToSelf");
              netMock->RequireReceive(&fixtureData.GetIoletId(), 1, 0, "IoletIdFromSelf");
              netMock->RequireSend(&expectedData.GetSiteType(), 1, 0, "SiteTypeToSelf");
              netMock->RequireReceive(&fixtureData.GetSiteType(), 1, 0, "SiteTypeFromSelf");

              netMock->RequireSend(exampleSite.GetWallDistances(),
                                   lb::lattices::D3Q15::NUMVECTORS - 1,
                                   0,
                                   "WallToSelf");
              netMock->RequireReceive(exampleSite.GetWallDistances(),
                                      lb::lattices::D3Q15::NUMVECTORS - 1,
                                      0,
                                      "WallFromSelf");
              netMock->RequireSend(&exampleSite.GetWallNormal(), 1, 0, "NormalToSelf");
              netMock->RequireReceive(&exampleSite.GetWallNormal(), 1, 0, "NormalFromSelf");
              manager->TransferNonFieldDependentInformation();
              netMock->ExpectationsAllCompleted();
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

              // begin by setting up mocks for the required site
              std::vector<int> countOfNeedsToZeroFromZero;
              countOfNeedsToZeroFromZero.push_back(1); // expectation
              std::vector<int> countOfNeedsFromZeroToZero;
              countOfNeedsFromZeroToZero.push_back(1); //fixture
              netMock->RequireSend(&countOfNeedsToZeroFromZero.front(), 1, 0, "CountToSelf");
              netMock->RequireReceive(&countOfNeedsFromZeroToZero.front(), 1, 0, "CountFromSelf");

              std::vector<site_t> needsShouldBeSentToSelf;
              std::vector<site_t> needsShouldBeReceivedFromSelf;
              needsShouldBeSentToSelf.push_back(targetGlobalOneDIdx); //expectation
              needsShouldBeReceivedFromSelf.push_back(targetGlobalOneDIdx); //fixture
              netMock->RequireSend(&needsShouldBeSentToSelf.front(), 1, 0, "NeedToSelf");
              netMock->RequireReceive(&needsShouldBeSentToSelf.front(), 1, 0, "NeedFromSelf");

              manager->RegisterNeededSite(targetGlobalOneDIdx);
              manager->ShareNeeds();
              netMock->ExpectationsAllCompleted();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);

              // It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

              netMock->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
                                   lb::lattices::D3Q15::NUMVECTORS,
                                   0,
                                   "IntersectionDataToSelf");

              std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
              netMock->RequireReceive(&(receivedFOld[0]),
                                      lb::lattices::D3Q15::NUMVECTORS,
                                      0,
                                      "IntersectionDataFromSelf");

              manager->TransferFieldDependentInformation();
              netMock->ExpectationsAllCompleted();

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

              // begin by setting up mocks for the required site
              std::vector<int> countOfNeedsToZeroFromZero;
              countOfNeedsToZeroFromZero.push_back(1); // expectation
              std::vector<int> countOfNeedsFromZeroToZero;
              countOfNeedsFromZeroToZero.push_back(1); //fixture
              netMock->RequireSend(&countOfNeedsToZeroFromZero.front(), 1, 0, "CountToSelf");
              netMock->RequireReceive(&countOfNeedsFromZeroToZero.front(), 1, 0, "CountFromSelf");

              std::vector<site_t> needsShouldBeSentToSelf;
              std::vector<site_t> needsShouldBeReceivedFromSelf;
              needsShouldBeSentToSelf.push_back(targetGlobalOneDIdx); //expectation
              needsShouldBeReceivedFromSelf.push_back(targetGlobalOneDIdx); //fixture
              netMock->RequireSend(&needsShouldBeSentToSelf.front(), 1, 0, "NeedToSelf");
              netMock->RequireReceive(&needsShouldBeSentToSelf.front(), 1, 0, "NeedFromSelf");

              manager->RegisterNeededSite(targetGlobalOneDIdx);
              manager->ShareNeeds();
              netMock->ExpectationsAllCompleted();

              // Now, transfer the data about that site.
              Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);
              // It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

              netMock->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
                                   lb::lattices::D3Q15::NUMVECTORS,
                                   0,
                                   "IntersectionDataToSelf");
              std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
              netMock->RequireReceive(&(receivedFOld[0]),
                                      lb::lattices::D3Q15::NUMVECTORS,
                                      0,
                                      "IntersectionDataFromSelf");

              net::phased::StepManager stepManager;
              stepManager.RegisterIteratedActorSteps(*manager, 0);
              net::phased::NetConcern netConcern = net::phased::NetConcern(*netMock);
              stepManager.RegisterCommsForAllPhases(netConcern);
              stepManager.CallActions();
              netMock->ExpectationsAllCompleted();

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
