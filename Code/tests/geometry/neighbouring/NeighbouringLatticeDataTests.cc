// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/neighbouring/NeighbouringDomain.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "lb/lattices/D3Q15.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/helpers/EqualitySiteData.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::geometry::neighbouring;

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "NeighbouringLatticeDataTests") {
      site_t dummyId = 54;
      auto&& data = latDat->GetNeighbouringData();
      auto&& dom = data.GetDomain();
      auto exampleSite = latDat->GetSite(24);

      SECTION("TestInsertAndRetrieveSiteData") {
	dom.GetSiteData(dummyId) = exampleSite.GetSiteData();
	REQUIRE(exampleSite.GetSiteData() == dom.GetSiteData(dummyId));
      }

      SECTION("TestInsertAndRetrieveDistance") {
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    dom.GetCutDistances(dummyId)[direction]=exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1);
	  }

	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    REQUIRE(exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1) == dom.GetCutDistance<lb::lattices::D3Q15>(dummyId, direction + 1));
	  }
      }
      
      SECTION("TestInsertAndRetrieveNormal") {
	dom.GetNormalToWall(dummyId) = exampleSite.GetWallNormal();
	REQUIRE(exampleSite.GetWallNormal() == dom.GetNormalToWall(dummyId));
      }

      SECTION("TestInsertAndRetrieveDistributions") {
	std::vector<distribn_t> distribution;
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
	  {
	    distribution.push_back(exampleSite.GetFOld<lb::lattices::D3Q15>()[direction]);
	  }

	data.GetDistribution(dummyId) = distribution;

	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
	  {
	    REQUIRE(exampleSite.GetFOld<lb::lattices::D3Q15>()[direction] ==
		    data.GetFOld(dummyId * lb::lattices::D3Q15::NUMVECTORS)[direction]);
	  }
      }

      SECTION("TestNeighbouringSite") {

	std::vector<distribn_t> distances;
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    distances.push_back(exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1));
	  }

	std::vector<distribn_t> distribution;
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
	  {
	    distribution.push_back(exampleSite.GetFOld<lb::lattices::D3Q15>()[direction]);
	  }
	data.SaveSite(dummyId,
		      distribution,
		      distances,
		      exampleSite.GetWallNormal(),
		      exampleSite.GetSiteData());

	auto neighbouringSite = data.GetSite(dummyId);
	REQUIRE(exampleSite.GetSiteData() == neighbouringSite.GetSiteData());
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    REQUIRE(exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1)
		    == neighbouringSite.GetWallDistance<lb::lattices::D3Q15>(direction + 1));
	  }
	REQUIRE(exampleSite.GetWallNormal() == neighbouringSite.GetWallNormal());
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS; direction++)
	  {
	    REQUIRE(exampleSite.GetFOld<lb::lattices::D3Q15>()[direction]
		    == neighbouringSite.GetFOld<lb::lattices::D3Q15>()[direction]);
	  }
      }
    }

  }
}

