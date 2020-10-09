
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <memory>

#include <catch2/catch.hpp>

#include "comm/MpiEnvironment.h"
#include "comm/AsyncConcern.h"
#include "geometry/LatticeData.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "timestep/TimeStepManager.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/helpers/MockCommunicator.h"
#include "tests/helpers/EqualitySiteData.h"

namespace hemelb
{
  namespace tests
  {
    
    using namespace hemelb::geometry::neighbouring;
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture, "NeighbouringDataManagerTests") {
      auto core_count = 1;
      auto current_core = 0;
      auto communicator = std::make_shared<helpers::MockCommunicator>(current_core, core_count);
      auto mockComms = std::dynamic_pointer_cast<helpers::MockCommunicator>(communicator);
      auto commQ = comm::Async::New(communicator);

      auto data = latDat->GetNeighbouringData();
      auto manager = NeighbouringDataManager{*latDat, data, commQ};
   

      SECTION("TestRegisterNeedOneProc") {
	// We imagine that unbeknownst to us, there is a site at
	// x=1,y=1,z=7 which for some unexplained reason we need to
	// know about.
	//
	// The client code has a duty to determine the global index
	// for the site, which it can do, using
	// LatticeData::GetGlobalNoncontiguousSiteIdFromGlobalCoords
	// and LatticeData::GetGlobalCoords. For our four cube test
	// case, we only have one proc, of course, so we pretend to
	// require a site from our own proc via this mechanism just to
	// verify it.

	util::Vector3D<site_t> exampleBlockLocalCoord(1, 2, 3);
	util::Vector3D<site_t> exampleBlockCoord(0, 0, 0);
	util::Vector3D<site_t> globalCoord = latDat->GetGlobalCoords(exampleBlockCoord,
								     exampleBlockLocalCoord);
	
	// to illustrate a typical use, we now add a displacement to
	// the global coords of the local site, which in a real
	// example would take it off-proc
	util::Vector3D<site_t> desiredGlobalCoord = globalCoord + util::Vector3D<site_t>{1,0,0};
	// desiredGC = (2,2,3) => id = (2*block + 2) * block + 3
	site_t desiredId = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(desiredGlobalCoord);
	REQUIRE(87 == desiredId);

	manager.RegisterNeededSite(desiredId);
	REQUIRE(1U == manager.GetNeededSites().size());
	REQUIRE(87 == manager.GetNeededSites()[0]);
      }

      SECTION("TestShareNeedsOneProc") {
	std::vector<site_t> neededsites = {43};
	manager.RegisterNeededSite(neededsites[0]);
             
	// The MapA2A needs a barrier but as the only needs are
	// local this is true immediately
	mockComms->RequireIbarrier([]() { return true; });

	// We then send this need to ourself
	mockComms->RequireSend(neededsites, 0, 0);
	mockComms->RequireRecv(neededsites, 0, 0);
      	manager.ShareNeeds();
      	mockComms->ExpectationsAllCompleted();

      	REQUIRE(manager.GetNeedsForProc(0).size() == 1U);
      	REQUIRE(manager.GetNeedsForProc(0).front() == 43);
      }

      using namespace geometry;
      SECTION("TestShareConstantDataOneProc") {
	std::vector<site_t> neededsites = {43};
	manager.RegisterNeededSite(neededsites[0]);

	mockComms->RequireIbarrier([]() { return true; });
	mockComms->RequireSend(neededsites, 0, 0);
	mockComms->RequireRecv(neededsites, 0, 0);
	manager.ShareNeeds();

	// Now, transfer the data about that site.
	Site < LatticeData > exampleSite
	  = latDat->GetSite(latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(43));
	// It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

	// We should send/receive the site data
	SiteData expectedData = exampleSite.GetSiteData();
	SiteData fixtureData = exampleSite.GetSiteData();

	mockComms->RequireSend(expectedData.GetWallIntersectionData(),
			       0, 0);
	mockComms->RequireRecv(fixtureData.GetWallIntersectionData(),
			       0, 0);
	mockComms->RequireSend(&expectedData.GetIoletIntersectionData(),
			       1,
			       0,
			       0);
	mockComms->RequireRecv(&fixtureData.GetIoletIntersectionData(),
			       1,
			       0,
			       0);
	mockComms->RequireSend(&expectedData.GetIoletId(), 1, 0, 0);
	mockComms->RequireRecv(&fixtureData.GetIoletId(), 1, 0, 0);
	mockComms->RequireSend(&expectedData.GetSiteType(), 1, 0, 0);
	mockComms->RequireRecv(&fixtureData.GetSiteType(), 1, 0, 0);

	mockComms->RequireSend(exampleSite.GetWallDistances(),
			       lb::lattices::D3Q15::NUMVECTORS - 1,
			       0,
			       0);
	mockComms->RequireRecv(exampleSite.GetWallDistances(),
			       lb::lattices::D3Q15::NUMVECTORS - 1,
			       0,
			       0);
	mockComms->RequireSend(&exampleSite.GetWallNormal(), 1, 0, 0);
	mockComms->RequireRecv(&exampleSite.GetWallNormal(), 1, 0, 0);
	manager.TransferNonFieldDependentInformation();
	      
	mockComms->ExpectationsAllCompleted();
      	NeighbouringSite transferredSite = data.GetSite(43);
      	REQUIRE(exampleSite.GetSiteData() == transferredSite.GetSiteData());
      	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++) {
      	  REQUIRE(exampleSite.GetWallDistances()[direction] == transferredSite.GetWallDistances()[direction]);
      	}
      	REQUIRE(exampleSite.GetWallNormal() == transferredSite.GetWallNormal());
      }

      SECTION("TestShareFieldDataOneProc") {
      	site_t targetGlobalOneDIdx = 43;
      	LatticeVector targetGlobalThreeDIdx = latDat->GetSiteCoordsFromSiteId(targetGlobalOneDIdx);
      	site_t targetLocalIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(targetGlobalOneDIdx);

      	for (unsigned int direction = 0; direction < 3; direction++)
      	  REQUIRE(1 == targetGlobalThreeDIdx[direction]);

	mockComms->RequireIbarrier([]() { return true; });
	std::vector<site_t> neededsites = {targetGlobalOneDIdx};
	mockComms->RequireSend(neededsites, 0, 0);
	mockComms->RequireRecv(neededsites, 0, 0);

	manager.RegisterNeededSite(targetGlobalOneDIdx);
	manager.ShareNeeds();

	// Now, transfer the data about that site.
	Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);

	// It should arrive in the NeighbouringDataManager, from the values sent from the localLatticeData

	mockComms->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
				 lb::lattices::D3Q15::NUMVECTORS,
				 0,
				 0);

	std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
	mockComms->RequireRecv(&(receivedFOld[0]),
				 lb::lattices::D3Q15::NUMVECTORS,
				 0,
				 0);
	
	manager.TransferFieldDependentInformation();
	mockComms->ExpectationsAllCompleted();

      	NeighbouringSite transferredSite = data.GetSite(targetGlobalOneDIdx);
      	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++) {
      	  REQUIRE(receivedFOld[direction] == transferredSite.GetFOld<lb::lattices::D3Q15> ()[direction]);
      	}
      }

      SECTION("TestShareFieldDataOneProcViaIterableAction") {
      	site_t targetGlobalOneDIdx = 43;
      	site_t targetLocalIdx = latDat->GetLocalContiguousIdFromGlobalNoncontiguousId(targetGlobalOneDIdx);

	mockComms->RequireIbarrier([]() { return true; });
	std::vector<site_t> neededsites = {targetGlobalOneDIdx};
	mockComms->RequireSend(neededsites, 0, 0);
	mockComms->RequireRecv(neededsites, 0, 0);

	manager.RegisterNeededSite(targetGlobalOneDIdx);
	manager.ShareNeeds();
	mockComms->ExpectationsAllCompleted();

	// Now, transfer the data about that site.
	Site < LatticeData > exampleSite = latDat->GetSite(targetLocalIdx);
	// It should arrive in the NeighbouringDataManager, from
	// the values sent from the localLatticeData

	mockComms->RequireSend(const_cast<distribn_t*> (exampleSite.GetFOld<lb::lattices::D3Q15> ()),
			       lb::lattices::D3Q15::NUMVECTORS,
			       0,
			       0);
	std::vector<distribn_t> receivedFOld(lb::lattices::D3Q15::NUMVECTORS, 53.0);
	mockComms->RequireRecv(&(receivedFOld[0]),
			       lb::lattices::D3Q15::NUMVECTORS,
			       0,
			       0);

	timestep::TimeStepManager stepManager(1);
	stepManager.AddToPhase(0, &manager);
	comm::AsyncConcern netConcern = comm::AsyncConcern(commQ);
	stepManager.AddToPhase(0, &netConcern);
	stepManager.DoStep();
	mockComms->ExpectationsAllCompleted();

      	NeighbouringSite transferredSite = data.GetSite(targetGlobalOneDIdx);
      	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++) {
      	  REQUIRE(receivedFOld[direction] == transferredSite.GetFOld<lb::lattices::D3Q15> ()[direction]);
      	}
      }

    }

  }
}

