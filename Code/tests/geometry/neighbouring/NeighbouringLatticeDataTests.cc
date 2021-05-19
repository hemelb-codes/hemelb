// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/neighbouring/NeighbouringLatticeData.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/helpers/EqualitySiteData.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::geometry::neighbouring;

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "NeighbouringLatticeDataTests") {
      site_t dummyId = 54;
      auto data = latDat->GetNeighbouringData();
      geometry::Site<geometry::LatticeData> exampleSite{latDat->GetSite(24)};

      SECTION("TestInsertAndRetrieveSiteData") {
	data.GetSiteData(dummyId) = exampleSite.GetSiteData();
	REQUIRE(exampleSite.GetSiteData() == data.GetSiteData(dummyId));
      }

      SECTION("TestInsertAndRetrieveDistance") {
	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    data.GetCutDistances(dummyId)[direction]=exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1);
	  }

	for (unsigned int direction = 0; direction < lb::lattices::D3Q15::NUMVECTORS - 1; direction++)
	  {
	    REQUIRE(exampleSite.GetWallDistance < lb::lattices::D3Q15 > (direction + 1) == data.GetCutDistance<lb::lattices::D3Q15>(dummyId, direction + 1));
	  }
      }
      
      SECTION("TestInsertAndRetrieveNormal") {
	data.GetNormalToWall(dummyId) = exampleSite.GetWallNormal();
	REQUIRE(exampleSite.GetWallNormal() == data.GetNormalToWall(dummyId));
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

	NeighbouringSite neighbouringSite = data.GetSite(dummyId);

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

