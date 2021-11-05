// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/CellEnergy.h"
#include "resources/Resource.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    static auto zero = Approx(0.0).margin(1e-8);

    // Helper type to do the tests on 
    template <typename FUNC, std::size_t... INDICES>
    struct EnergyChecker {
      MeshData const& mesh;
      FUNC func;

      template <std::size_t... T>
      EnergyChecker(MeshData const& m, FUNC f) : mesh{m}, func{f} {
      }

      void operator()(LatticePosition const &normal, size_t node) {
	std::vector<LatticeForceVector> forces(4, LatticeForceVector(0, 0, 0));
	double const epsilon(1e-5);
	MeshData newmesh(mesh);
	newmesh.vertices[node] += normal * epsilon;
	// Here we call the function = have to push in the indices to
	// use for the facetBending.
	LatticeEnergy const deltaE(func(newmesh.vertices, mesh, INDICES..., 1.0, forces));
	REQUIRE(deltaE / epsilon == zero);
	for (auto& f: forces)
	  REQUIRE(f == ApproxV(LatticeForceVector::Zero()));
	
      }
    };

    // Helper so the compiler will deduce FUNC for us.
    template <std::size_t... INDICES, typename FUNC>
    EnergyChecker<FUNC, INDICES...> MakeEnergyChecker(MeshData const& m, FUNC f) {
      return EnergyChecker<FUNC, INDICES...>{m, f};
    }

    // Quick mildly hacky macro to "pass an overload set" as the function argument above
#define LIFT(func) [](auto... x){ return func(x...);}

    // Tests certain node movement that do not result in forces/energy
    TEST_CASE_METHOD(BasisFixture, "GradientKernTests", "[redblood]") {
      SECTION("testBending") {
	auto noBending = MakeEnergyChecker<0, 3>(mesh, LIFT(facetBending));
	// Single nodes
	noBending(LatticeForceVector(1, 0, 0).GetNormalised(), 0);
	noBending(LatticeForceVector(0, 1, 0).GetNormalised(), 0);
	noBending(LatticeForceVector(1, 1, -2).GetNormalised(), 3);
	noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 3);
	// Common nodes
	noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 1);
	noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 2);
      }

      SECTION("testVolume") {
	auto noVolume = MakeEnergyChecker<>(mesh, LIFT(volumeEnergy));
	noVolume(LatticePosition(-1, 1, 0).GetNormalised(), 0);
	noVolume(LatticePosition(1, 1, -2).GetNormalised(), 0);
	noVolume(LatticePosition(0, 1, 0).GetNormalised(), 1);
	noVolume(LatticePosition(0, 0, 1).GetNormalised(), 1);
	noVolume(LatticePosition(1, 0, 0).GetNormalised(), 2);
	noVolume(LatticePosition(0, 0, 1).GetNormalised(), 2);
	noVolume(LatticePosition(1, 0, 0).GetNormalised(), 3);
	noVolume(LatticePosition(0, 1, 0).GetNormalised(), 3);
      }

      SECTION("testSurface") {
	// Get rid of some surfaces first, to simplify kernel
	MeshData const save(mesh);
	mesh.facets.resize(2);
	mesh.facets.back() = save.facets.back();
	auto noSurface = MakeEnergyChecker<>(mesh, LIFT(surfaceEnergy));
	noSurface(LatticePosition(0, 0, 1).GetNormalised(), 0);
	noSurface(LatticePosition(1, 1, 1).GetNormalised(), 3);
	noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 0);
	noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 3);
      }
    }
  }
}
