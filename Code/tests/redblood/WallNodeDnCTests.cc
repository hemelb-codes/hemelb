// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/WallCellPairIterator.h"
#include "lb/lattices/D3Q15.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "WallNodeDnCTests", "[redblood]") {
	
      LatticeDistance const cutoff = 3.0;
      LatticeDistance const halo = 0.5;
      using Lattice = lb::lattices::D3Q15;

      std::unique_ptr<tests::FourCubeLatticeData> latticeData(FourCubeLatticeData::Create(Comms(), 27 + 2));

      for (site_t i(0); i < latticeData->GetLocalFluidSiteCount(); ++i) {
	auto const site = latticeData->GetSite(i);
	if (not site.IsWall()) {
	  continue;
	}
	for (Direction d(0); d < Lattice::NUMVECTORS; ++d) {
	  if (site.HasWall(d)) {
	    latticeData->SetBoundaryDistance(i, d, 0.5);
	  }
	}
      }

      SECTION("testWallNodeDnC") {
	auto const dnc = createWallNodeDnC<Lattice>(*latticeData, cutoff, halo);

	// Checking the middle of the world first: No wall nodes
	auto const center = dnc.equal_range(LatticePosition(13, 13, 13));
	REQUIRE(center.first == center.second);

	// Checking that upper side has 3*3*5 wall sites,
	// 3*3 because that's the size of a DnC box
	// 5 because there are five possible directions in D3Q15 pointing towards the wall
	auto upper = dnc.equal_range(LatticePosition(0, 3.5 * 3, 3.5 * 3));
	REQUIRE(45l == std::distance(upper.first, upper.second));
	for (; upper.first != upper.second; ++upper.first) {
	  REQUIRE(upper.first->second.nearBorder bitand size_t(Borders::CENTER));
	  bool cond1 = (upper.first->second.nearBorder bitand size_t(Borders::BOTTOM)) != 0;
	  bool cond2 = upper.first->second.node.x < halo;
	  REQUIRE(cond1 == cond2);
	}
      }
    }

  }
}
