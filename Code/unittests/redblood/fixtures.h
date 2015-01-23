//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H
#define HEMELB_UNITTESTS_REDBLOOD_FIXTURES_H

#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"

namespace hemelb { namespace unittests {

class BasisFixture : public CppUnit::TestFixture {
  public:
    void setUp() {
      // facets at something degrees from one another
      mesh.vertices.push_back(LatticePosition(0, 0, 0));
      mesh.vertices.push_back(LatticePosition(1, 0, 0));
      mesh.vertices.push_back(LatticePosition(0, 1, 0));
      mesh.vertices.push_back(LatticePosition(0, 0, 1));


      redblood::MeshData::Facet indices;
      indices[0] = 0; indices[1] = 1; indices[2] = 2;
      mesh.facets.push_back(indices);
      indices[0] = 0; indices[1] = 2; indices[2] = 3;
      mesh.facets.push_back(indices);
      indices[0] = 0; indices[1] = 3; indices[2] = 1;
      mesh.facets.push_back(indices);
      indices[0] = 1; indices[1] = 3; indices[2] = 2;
      mesh.facets.push_back(indices);
    }

  protected:
    hemelb::redblood::MeshData mesh;
};

class EnergyVsGradientFixture : public BasisFixture {

  public:
  template<class ENERGY, class GRADIENT>
    void energyVsForces(
        ENERGY const &_energy,
        GRADIENT const &_gradient,
        LatticePosition const &_dir,
        size_t _node, double _epsilon = 1e-8) {

      std::vector<LatticeForceVector> forces(4, LatticeForceVector(0,0,0));
      PhysicalEnergy const firstE(_gradient(mesh, forces));

      redblood::MeshData newmesh(mesh);
      newmesh.vertices[_node] += _dir * _epsilon;
      PhysicalEnergy const deltaE(_energy(newmesh) - firstE);

      double const tolerance(
          std::max(std::abs((deltaE / _epsilon) * 1e-4), 1e-8)
      );
      CPPUNIT_ASSERT(
          helpers::is_zero(
            forces[_node].Dot(_dir) + (deltaE / _epsilon),
            tolerance
          )
      );
    }

  template<class ENERGY, class GRADIENT>
    void energyVsForces(ENERGY const &_energy,
        GRADIENT const &_gradient,
        double _epsilon = 1e-8) {

      for(size_t node(0); node < mesh.vertices.size(); ++node) 
        for(size_t i(0); i < 3; ++i)
          energyVsForces(
            _energy, _gradient,
            LatticePosition(i==0, i==1, i==2),
            node, _epsilon
          );
    }

  template<class BOTH>
    void energyVsForces(BOTH const &_both, double _epsilon = 1e-8) {
      energyVsForces(_both, _both, _epsilon);
    }
};

class SquareDuctTetrahedronFixture : public helpers::FourCubeBasedTestFixture {
  public:
    SquareDuctTetrahedronFixture()
      : FourCubeBasedTestFixture(), mesh(redblood::tetrahedron()) {}

    void setUp() {
      FourCubeBasedTestFixture::setUp();
      mesh = redblood::refine(initial_mesh(), refinement());
      mesh *= Dimensionless(CubeSize() - 3) * 0.5;
      mesh += Dimensionless(CubeSize()) * 0.5;
    }

  protected:
    // Functions to parameterize the fizture
    virtual size_t CubeSize() const { return 32 + 2; }
    virtual size_t refinement() const { return 3; }
    virtual redblood::Mesh initial_mesh() const {
      return redblood::tetrahedron();
    }

    hemelb::redblood::Cell mesh;
};

template<class CELLTYPE=redblood::Cell>
  redblood::CellContainer TwoPancakeSamosas(PhysicalDistance _cutoff) {
    redblood::CellContainer cells;
    redblood::Mesh pancake = redblood::pancakeSamosa();
    pancake += LatticePosition(1, 1, 1) * _cutoff * 0.5;
    // safer to clone so cells has its own copy
    cells.emplace_back(std::make_shared<CELLTYPE>(pancake.clone()));
    pancake += LatticePosition(3, 0, 1) * _cutoff;
    cells.emplace_back(std::make_shared<CELLTYPE>(pancake.clone()));

    return cells;
  }
}}
#endif
