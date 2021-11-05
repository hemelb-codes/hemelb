// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "resources/Resource.h"
#include "redblood/Facet.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    // Tests functionality that is *not* part of the HemeLB API
    // Checks that we know how to compute geometric properties between facets
    // However, HemeLB only cares about energy and forces
    TEST_CASE_METHOD(BasisFixture, "FacetTests", "[redblood]") {
      auto main = std::make_shared<Facet>(mesh, 0);
      auto neighbor = std::make_shared<Facet>(mesh, 3);

      SECTION("testNormal") {
	{
	  LatticePosition const actual = main->normal();
	  LatticePosition const expected(0, 0, -1);
	  REQUIRE(actual == ApproxV(expected));
	}
	{
	  LatticePosition const actual = neighbor->normal();
	  LatticePosition const expected(1, 1, 1);
	  REQUIRE(actual == ApproxV(expected));
	}
      }

      SECTION("testUnitNormal") {
	{
	  LatticePosition const actual = main->unitNormal();
	  LatticePosition const expected(0, 0, -1);
	  REQUIRE(actual == ApproxV(expected));
	}
	{
	  LatticePosition const actual = neighbor->unitNormal();
	  LatticePosition const expected(1, 1, 1);
	  REQUIRE(actual == ApproxV(expected.GetNormalised()));
	}
      }

      SECTION("testAngle") {
	Angle const actual0 = angle(*main, *neighbor);
	REQUIRE(actual0 == Approx(std::acos(-1.0 / std::sqrt(3.))));

	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	Angle const actual1 = angle(*main, *neighbor);
	REQUIRE(actual1 == Approx(3.0 * PI / 4.0));
	mesh.vertices.back()[1] = 1e0;
      }

      SECTION("testCommonNodes") {
	IndexPair const nodes = commonNodes(*main, *neighbor);
	REQUIRE(nodes.first == 1);
	REQUIRE(nodes.second == 2);
	LatticePosition const edge = (*main)(nodes.first, nodes.second);
	REQUIRE(edge == ApproxV(mesh.vertices[1] - mesh.vertices[2]));
      }
      SECTION("testSingleNodes") {
	IndexPair const nodes = singleNodes(*main, *neighbor);
	// indices into the list of vertices of each facet,
	// not into list of all vertices.
	REQUIRE(nodes.first == 0);
	REQUIRE(nodes.second == 1);
      }

      SECTION("testOrientedAngle") {
	Angle const actual0 = orientedAngle(*main, *neighbor);
	REQUIRE(Approx(std::acos(-1.0 / std::sqrt(3.))).margin(1e-8) == actual0);

	// simpler angle
	mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
	Angle const actual1 = orientedAngle(*main, *neighbor);
	REQUIRE(Approx(3.0 * PI / 4.0).margin(1e-8) == actual1);

	// change orientation <==> negative angle
	mesh.facets.front()[1] = 2;
	mesh.facets.front()[2] = 1;
	Angle const actual2 = orientedAngle(Facet(mesh, 0), Facet(mesh, 3));
	REQUIRE(Approx(-PI / 4.0).margin(1e-8) == actual2);
	mesh.facets.front()[1] = 1;
	mesh.facets.front()[2] = 2;
	mesh.vertices.back()[1] = 1e0;
      }


    }
  }
}
