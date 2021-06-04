// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "redblood/CellEnergy.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    TEST_CASE_METHOD(BasisFixture, "EnergyTests", "[redblood]") {
      MeshData original = mesh;
      std::vector<LatticeForceVector> forces{4, LatticeForceVector::Zero()};

      SECTION("testBending") {
	// No difference between original and current mesh
	// Hence energy is zero
	LatticeEnergy const actual0(facetBending(mesh.vertices, original, 0, 3, 1e0));
	REQUIRE(Approx(0e0) == actual0);

	// Now modify mesh and check "energy" is square of angle difference
	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	LatticeEnergy const actual1(facetBending(mesh.vertices, original, 0, 3, 1e0));
	mesh.vertices.back()[2] = 1e0;

	LatticeEnergy const expected(std::pow( (PI / 4e0 - std::acos(1. / std::sqrt(3.))), 2));
	REQUIRE(Approx(std::sqrt(3.) * expected) == actual1);
      }

      SECTION("testVolume") {
	// No difference between original and current mesh
	// Hence energy is zero
	LatticeEnergy const actual0(volumeEnergy(mesh.vertices, original, 1e0));
	REQUIRE(Approx(0.0) == actual0);

	// Now modify mesh and check "energy" is square of volume diff
	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	LatticeEnergy const actual1(volumeEnergy(mesh.vertices,
						 original,
						 2.0 * volume(original)));

	LatticeEnergy const deltaV(volume(mesh) - volume(original));
	REQUIRE(actual1 == Approx(deltaV * deltaV));
	mesh.vertices.back()[2] = 1e0;
      }

      SECTION("testSurface") {
	// No difference between original and current mesh
	// Hence energy is zero
	LatticeEnergy const actual0(surfaceEnergy(mesh.vertices, original, 1e0));
	REQUIRE(Approx(0) == actual0);

	// Now modify mesh and check "energy" is square of volume diff
	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	LatticeEnergy const actual1(surfaceEnergy(mesh.vertices,
						  original,
						  2.0 * area(original)));

	LatticeEnergy const deltaS(area(mesh) - area(original));
	REQUIRE(actual1 == Approx(deltaS * deltaS));
	mesh.vertices.back()[2] = 1e0;
      }

      SECTION("testStrain") {
	// No difference between original and current mesh
	// Hence energy is zero
	LatticeEnergy const actual0(strainEnergy(mesh.vertices, original, 1e0, 2e0));
	REQUIRE(actual0 == Approx(0.0));

	// Now modify mesh and check "energy" is square of volume diff
	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	LatticeEnergy const actual1(strainEnergy(mesh.vertices, original, 1e0, 2e0));

	LatticeEnergy const regression(0.0865562612162);
	REQUIRE(actual1 == Approx(regression));
	mesh.vertices.back()[2] = 1e0;
      }


    }
  }
}
