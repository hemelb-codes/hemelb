// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REDBLOOD_FIXTURES_H
#define HEMELB_TESTS_REDBLOOD_FIXTURES_H

#include <memory>
#include <catch2/catch.hpp>

#include "redblood/types_fwd.h"
#include "redblood/Cell.h"
#include "redblood/FlowExtension.h"
#include "redblood/Mesh.h"

#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb {

  namespace tests {
    template<class CELLTYPE = redblood::Cell>
    redblood::CellContainer TwoPancakeSamosas(LatticeDistance cutoff)
    {
      redblood::CellContainer cells;
      redblood::Mesh pancake = redblood::pancakeSamosa();
      pancake += LatticePosition(1, 1, 1) * cutoff * 0.5;
      // safer to clone so cells has its own copy
      cells.insert(std::make_shared<CELLTYPE>(pancake.clone()));
      pancake += LatticePosition(3, 0, 1) * cutoff;
      cells.insert(std::make_shared<CELLTYPE>(pancake.clone()));

      return cells;
    }

    class BasisFixture {
    public:
      BasisFixture();
    protected:
      hemelb::redblood::MeshData mesh;
    };

    class EnergyVsGradientFixture : public BasisFixture
    {
    public:
      template<class ENERGY, class GRADIENT>
      void energyVsForces(ENERGY const &energy, GRADIENT const &gradient,
			  LatticePosition const &dir, size_t node, double epsilon = 1e-8)
      {
	std::vector<LatticeForceVector> forces(4, LatticeForceVector(0, 0, 0));
	LatticeEnergy const firstE(gradient(mesh, forces));

	redblood::MeshData newmesh(mesh);
	newmesh.vertices[node] += dir * epsilon;
	LatticeEnergy const deltaE(energy(newmesh) - firstE);

	double const tolerance(std::max(std::abs( (deltaE / epsilon) * 1e-4), 1e-8));
	REQUIRE(Approx(- (deltaE / epsilon)).margin(tolerance) == forces[node].Dot(dir));
      }

      template<class ENERGY, class GRADIENT>
      void energyVsForces(ENERGY const &energy, GRADIENT const &gradient, double epsilon = 1e-8)
      {
	for (size_t node(0); node < mesh.vertices.size(); ++node)
	  for (size_t i(0); i < 3; ++i)
	    energyVsForces(energy,
			   gradient,
			   LatticePosition(i == 0, i == 1, i == 2),
			   node,
			   epsilon);
      }

      template<class BOTH>
      void energyVsForces(BOTH const &both, double epsilon = 1e-8)
      {
	energyVsForces(both, both, epsilon);
      }
    };

    class SquareDuctTetrahedronFixture : public helpers::FourCubeBasedTestFixture<32>
    {
    protected:
      size_t refinement = 3;
      hemelb::redblood::Cell mesh;

    public:
      SquareDuctTetrahedronFixture(redblood::Mesh const & initial_mesh, size_t refinement);
    };

    class FlowExtensionFixture
    {
    public:
      FlowExtensionFixture();
    protected:
      hemelb::redblood::FlowExtension flowExt;
    };

  }
}

#endif
