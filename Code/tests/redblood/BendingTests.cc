// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iomanip>
#include <catch2/catch.hpp>

#include "redblood/Cell.h"
#include "redblood/CellEnergy.h"
#include "redblood/Facet.h"
#include "redblood/Mesh.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    TEST_CASE("BendingTests") {
      auto vertices = MeshData::Vertices{ { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
      MeshData mesh;
      mesh.vertices = vertices;
      mesh.facets = MeshData::Facets{ { { 0, 1, 2 } }, { { 1, 3, 2 } } };
      std::vector<LatticeForceVector> forces{4, LatticeForceVector(0, 0, 0)};
      LatticeForce const moduli = 1e0;

      auto z10 = Approx(0.0).margin(1e-10);
      auto z8 = Approx(0.0).margin(1e-8);

      auto energyAndForces = [&] () -> LatticeEnergy {
	return facetBending(vertices, mesh, 0, 1, moduli, forces);
      };
      auto energy = [&]() -> LatticeEnergy {
	return facetBending(vertices, mesh, 0, 1, moduli);
      };

      auto numericalForcesSingle = [&] (LatticePosition const &direction, int node) {
      	auto const theta0 = orientedAngle(Facet(mesh, 0), Facet(mesh, 1));
      	auto const theta = orientedAngle(Facet(vertices, mesh.facets, 0),
      					 Facet(vertices, mesh.facets, 1));
      	auto const epsilon = 1e-4;
      	std::fill(forces.begin(), forces.end(), LatticeForceVector::Zero());
      	auto const oldPos = vertices[node];

      	auto const e0 = energyAndForces();
      	if (std::abs(theta - theta0) < 1e-8) {
	  REQUIRE(z8 == forces[node].Dot(direction));
	  return;
	}

	vertices[node] += direction * epsilon;
	auto const e1 = energy();
	auto const tol = std::max(std::abs( (e1 - e0) / epsilon) * 1e-2, 1e-5);
	REQUIRE(Approx(- (e1 - e0) / epsilon).margin(tol) == forces[node].Dot(direction));
	vertices[node] = oldPos;
      };

      auto numericalForces = [&]() {
	for (int j(0); j < 6; ++j) {
	  LatticePosition const direction(j == 0 ?
					  1 :
					  (j == 1 ?
					   -1 :
					   (j >= 6 ?
					    random() - 0.5 :
					    0)),
					  j == 2 ?
					  1 :
					  (j == 3 ?
					   -1 :
					   (j >= 6 ?
					    random() - 0.5 :
					    0)),
					  j == 4 ?
					  1 :
					  (j == 5 ?
					   -1 :
					   (j >= 6 ?
					    random() - 0.5 :
					    0)));
	  numericalForcesSingle(direction.GetNormalised(), 0);
	  numericalForcesSingle(direction.GetNormalised(), 1);
	  numericalForcesSingle(direction.GetNormalised(), 2);
	  numericalForcesSingle(direction.GetNormalised(), 3);
	}
      };

          auto bending = [&](LatticePosition const &v, Dimensionless theta) -> LatticePosition
          {
            auto const &a0 = mesh.vertices[1];
            auto const &a1 = mesh.vertices[2];
            return rotationMatrix(a1 - a0, theta) * (v - a0) + a0;
          };
          // LatticePosition twisting(LatticePosition const &v, Dimensionless theta) const
          // {
          //   auto const &a0 = mesh.vertices[2];
          //   return rotationMatrix(LatticePosition(1, 1, 1), theta) * (v - a0) + a0;
          // }

      SECTION("testNoBendingNoNothing") {
	REQUIRE(energyAndForces() == z10);
	for (auto const& force : forces) {
	  REQUIRE(force.x == z10);
	  REQUIRE(force.y == z10);
	  REQUIRE(force.z == z10);
	}
      }

      SECTION("testEnergy") {
	for (auto const theta : { 1e-2, 2e-2, 3e-2 }) {
	  std::fill(forces.begin(), forces.end(), LatticeForceVector::Zero());
	  vertices.back() = bending(mesh.vertices.back(), theta);
	  REQUIRE(z10(std::sqrt(3.) * theta * theta * moduli) == energy());
	}
      }

      SECTION("testStuff") {
	// modify template geometry to be almost planar
	auto const theta = 1e-2;
	mesh.vertices.back() = bending(LatticePosition(1, 1, 0), theta);
	// Now go to flat geometry, check energy and forces
	vertices.back() = LatticePosition(1, 1, 0);
	REQUIRE(z10(std::sqrt(3.) * moduli * theta * theta) == energyAndForces());
	for (auto const &force : forces) {
	  REQUIRE(z8 == force.x);
	  REQUIRE(z8 == force.y);
	  REQUIRE(z8 == force.z);
	}
      }

      SECTION("testNumericalForces") {
	numericalForces();
      }

      SECTION("testInflate") {
	for (auto &vertex : vertices)
	  {
	    vertex = vertex * 1.5;
	  }
	numericalForces();
      }

      SECTION("testConvexities") {

	auto with_angles = [&](double theta0, double theta)
	  {
	    mesh.vertices.back() = bending(LatticePosition(1, 1, 0), theta0);
	    vertices.back() = bending(LatticePosition(1, 1, 0), theta);
	    numericalForces();
	  };
	// Explores all sorts of configurations with different concavities, etc
	for (auto const theta0 : { 2e-1, -2e-1 })
	  {
	    for (auto const theta : { 1e-1, -1e-1, 3e-1, -3e-1 })
              {
                if (std::abs(theta0 - theta) > 1e-8)
		  {
		    with_angles(theta0, theta);
		  }
              }
	  }
      }

      SECTION("testSwapFacets") {
	std::swap(mesh.facets[0], mesh.facets[1]);
	numericalForces();
      }

      SECTION("testRotateNodeInFacet") {
	std::swap(mesh.facets[0][0], mesh.facets[0][1]);
	std::swap(mesh.facets[0][0], mesh.facets[0][2]);
	numericalForces();

	std::swap(mesh.facets[1][0], mesh.facets[1][1]);
	std::swap(mesh.facets[1][0], mesh.facets[1][2]);
	numericalForces();
      }

      SECTION("testMoveNodes") {
	for (auto &vertex : vertices) {
	  for (int j(0); j < 6; ++j) {
	    auto const epsilon = 1e-1;
	    LatticePosition const direction(j == 0 ?
					    1 :
					    (j == 1 ?
					     -1 :
					     0),
					    j == 2 ?
					    1 :
					    (j == 3 ?
					     -1 :
					     0),
					    j == 4 ?
					    1 :
					    (j == 5 ?
					     -1 :
					     0));
	    auto const old = vertex;
	    vertex += direction * epsilon;
	    numericalForces();
	    vertex = old;
	  }
	}
      }

    }
  }
}
