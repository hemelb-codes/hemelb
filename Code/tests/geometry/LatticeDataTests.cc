// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "geometry/LatticeData.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::geometry;
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "LatticeDataTests") {
      SECTION("TestConvertGlobalId") {
	// Not really a very good test to use a one-proc geometry.  We
	// need to create a sixteen-cube lattice data test fixture to
	// test this kind of thing properly, but we can at least smoke
	// test things with four cube.
	util::Vector3D<site_t> exampleCoord(1, 2, 3);
	site_t id = latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(exampleCoord);
	site_t blockSize = latDat->GetBlockSize();
	REQUIRE(site_t( (1 * blockSize + 2) * blockSize + 3) == id);
	util::Vector3D<site_t> resultCoord;
	latDat->GetGlobalCoordsFromGlobalNoncontiguousSiteId(id, resultCoord);
	REQUIRE(resultCoord == exampleCoord);
      }

      SECTION("TestGetProcFromGlobalId") {
	// Again, we need to work up a way to mock a multi-processor
	// situation to test this properly.
	REQUIRE(latDat->ProcProvidingSiteByGlobalNoncontiguousId(43) == 0);
      }
    }
  }
}

