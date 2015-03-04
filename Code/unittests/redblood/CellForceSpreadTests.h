//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELL_FORCE_SPREAD_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELL_FORCE_SPREAD_TESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/GridAndCell.h"
#include "redblood/facet.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellForceSpreadTests : public SquareDuctTetrahedronFixture
      {
          CPPUNIT_TEST_SUITE (CellForceSpreadTests);
          CPPUNIT_TEST (testIsZeroFarFromMembrane);
          CPPUNIT_TEST (testIsSymmetric);
          CPPUNIT_TEST (testIsIncreasing);
          CPPUNIT_TEST (testIsLinear);CPPUNIT_TEST_SUITE_END();

          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> Kernel;

        public:
          void setUp();
          void setUpForces();

          // Moves fixture so barycenter is at given position and computes spread
          // force at origin
          LatticeForceVector force_at_center(LatticePosition const &position);
          void testIsZeroFarFromMembrane();
          void testIsSymmetric();
          void testIsIncreasing();
          void testIsLinear();

        protected:
          std::vector<LatticeForceVector> forces;
          LatticePosition direction, intensity, center;
          virtual Mesh initial_mesh() const
          {
            return triangleMesh();
          }
          virtual size_t refinement() const
          {
            return 3;
          }
      };

      LatticeForceVector CellForceSpreadTests::force_at_center(LatticePosition const &position)
      {
        mesh += position - mesh.GetBarycenter();
        helpers::ZeroOutForces(latDat);
        details::spreadForce2Grid(
            std::shared_ptr<CellBase>(&mesh, [](CellBase*){}),
            details::SpreadForces(forces, *latDat),
            stencil::types::FOUR_POINT);
        return latDat->GetSite(center).GetForce();
      }

      void CellForceSpreadTests::setUp()
      {
        center = LatticePosition(CubeSize() / 2, CubeSize() / 2, CubeSize() / 2);
        SquareDuctTetrahedronFixture::setUp();
        setUpForces();
      }
      void CellForceSpreadTests::setUpForces()
      {
        direction = LatticePosition(1, 2, 3);
        intensity = LatticePosition(3, 2, 1).Normalise();
        forces.resize(mesh.GetNumberOfNodes());
        typedef MeshData::Vertices::const_iterator const_iterator;
        const_iterator i_vertex = mesh.GetVertices().begin();
        const_iterator const i_end = mesh.GetVertices().end();
        std::vector<LatticeForceVector>::iterator i_force = forces.begin();

        for (; i_vertex != i_end; ++i_vertex, ++i_force)
        {
          *i_force = direction * i_vertex->Dot(intensity);
        }
      }

      void CellForceSpreadTests::testIsZeroFarFromMembrane()
      {
        // Very far away from pancake samosa
        LatticeForceVector const faraway = force_at_center(center + LatticePosition(0, 0, 3));
        CPPUNIT_ASSERT(helpers::is_zero(faraway));

        LatticeForceVector const justToFar = force_at_center(center + LatticePosition(0, 0, 2));
        CPPUNIT_ASSERT(helpers::is_zero(justToFar));

        LatticeForceVector const justInside = force_at_center(center
            + LatticePosition(0, 0, 2.0 - 1e-4));
        CPPUNIT_ASSERT(not helpers::is_zero(justInside));
        CPPUNIT_ASSERT(justInside[0] > 1e-8 and justInside[0] < 1e-4);
      }

      void CellForceSpreadTests::testIsSymmetric()
      {
        size_t const N(10);

        for (size_t i(0); i < N; ++i)
        {
          LatticePosition const displacement(0, 0, 2e0 - Dimensionless(i * 2) / Dimensionless(N));
          LatticeForceVector const left = force_at_center(center + displacement);
          LatticeForceVector const right = force_at_center(center + displacement);
          CPPUNIT_ASSERT(helpers::is_zero(left - right));
        }
      }

      void CellForceSpreadTests::testIsIncreasing()
      {
        size_t const N(10);
        LatticeForceVector last(0, 0, 0);

        for (size_t i(0); i < N; ++i)
        {
          LatticePosition const displacement(0,
                                             0,
                                             2e0 - Dimensionless(i * 2 + 1) / Dimensionless(N));
          LatticeForceVector const current = force_at_center(center + displacement);
          CPPUNIT_ASSERT(current[0] > last[0]);
          CPPUNIT_ASSERT(current[1] > last[1]);
          CPPUNIT_ASSERT(current[2] > last[2]);
          last = current;
        }
      }

      void CellForceSpreadTests::testIsLinear()
      {
        size_t const N(5);
        mesh = Cell(refine(MeshData { mesh.GetVertices(), mesh.GetFacets() }, 4));
        setUpForces();
        // x0, x1 should be further than 2 from the edges
        // Only linear if samosa appears as infinite plane
        // with sufficiently dense vertices
        LatticePosition const x0(center[0], center[1] - 0.5, center[2] - 0.1), x1(center[0],
                                                                                  center[1] + 0.5,
                                                                                  center[2] - 0.1);
        LatticeForceVector const v0(force_at_center(x0)), v1(force_at_center(x1));

        LatticeForceVector const a( (v1 - v0) / (direction.Dot(x1) - direction.Dot(x0)));
        Dimensionless const tolerance(std::max(std::max(std::abs( (v0 - v1)[0]),
                                                        std::abs( (v0 - v1)[1])),
                                               std::abs( (v0 - v1)[2])) * 1e-3);

        for (size_t i(0); i < N; ++i)
        {
          LatticePosition const x = (x1 - x0) * (Dimensionless(i + 1) / Dimensionless(N + 2)) + x0;
          LatticeForceVector const expected(a * (direction.Dot(x) - direction.Dot(x0)) + v0);
          CPPUNIT_ASSERT(helpers::is_zero(expected - force_at_center(x), tolerance));
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellForceSpreadTests);
    }
  }
}

#endif  // ONCE
