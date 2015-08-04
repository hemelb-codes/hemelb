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
#include "redblood/FlowExtension.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace unittests
  {
    class BasisFixture : public CppUnit::TestFixture
    {
      public:
        void setUp()
        {
          // facets at something degrees from one another
          mesh.vertices.push_back(LatticePosition(0, 0, 0));
          mesh.vertices.push_back(LatticePosition(1, 0, 0));
          mesh.vertices.push_back(LatticePosition(0, 1, 0));
          mesh.vertices.push_back(LatticePosition(0, 0, 1));

          redblood::MeshData::Facet indices;
          indices[0] = 0;
          indices[1] = 1;
          indices[2] = 2;
          mesh.facets.push_back(indices);
          indices[0] = 0;
          indices[1] = 2;
          indices[2] = 3;
          mesh.facets.push_back(indices);
          indices[0] = 0;
          indices[1] = 3;
          indices[2] = 1;
          mesh.facets.push_back(indices);
          indices[0] = 1;
          indices[1] = 3;
          indices[2] = 2;
          mesh.facets.push_back(indices);
          redblood::orientFacets(mesh);
        }

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
          PhysicalEnergy const firstE(gradient(mesh, forces));

          redblood::MeshData newmesh(mesh);
          newmesh.vertices[node] += dir * epsilon;
          PhysicalEnergy const deltaE(energy(newmesh) - firstE);

          double const tolerance(std::max(std::abs( (deltaE / epsilon) * 1e-4), 1e-8));
          CPPUNIT_ASSERT_DOUBLES_EQUAL(-(deltaE / epsilon), forces[node].Dot(dir), tolerance);
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

    class SquareDuctTetrahedronFixture : public helpers::FourCubeBasedTestFixture
    {
      public:
        SquareDuctTetrahedronFixture() :
            FourCubeBasedTestFixture(), mesh(redblood::tetrahedron())
        {
        }

        void setUp()
        {
          FourCubeBasedTestFixture::setUp();
          mesh = redblood::refine(initial_mesh(), refinement());
          mesh *= Dimensionless(CubeSize() - 3) * 0.5;
          mesh += Dimensionless(CubeSize()) * 0.5;
        }

      protected:
        // Functions to parameterize the fizture
        virtual size_t CubeSize() const
        {
          return 32 + 2;
        }
        virtual size_t refinement() const
        {
          return 3;
        }
        virtual redblood::Mesh initial_mesh() const
        {
          return redblood::tetrahedron();
        }

        hemelb::redblood::Cell mesh;
    };

    class FlowExtensionFixture : public CppUnit::TestFixture
    {
      public:
        void setUp()
        {
          flowExt.normal = util::Vector3D<LatticeDistance>(1.0, 0.0, 0.0);
          flowExt.origin = util::Vector3D<LatticeDistance>(0.0, 0.0, 0.0);
          flowExt.length = 10.0;
          flowExt.radius = 1.0;
        }
      protected:
        hemelb::redblood::FlowExtension flowExt;
    };

    template<class CELLTYPE = redblood::Cell>
    redblood::CellContainer TwoPancakeSamosas(PhysicalDistance cutoff)
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
  }
}
#endif
