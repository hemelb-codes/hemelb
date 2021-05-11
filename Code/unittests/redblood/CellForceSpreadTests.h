// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLFORCESPREADTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLFORCESPREADTESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/GridAndCell.h"
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
          CPPUNIT_TEST (testIsZeroFarFromMembrane<stencil::FourPoint> );
          CPPUNIT_TEST (testIsZeroFarFromMembrane<stencil::CosineApprox> );
          CPPUNIT_TEST (testIsZeroFarFromMembrane<stencil::ThreePoint> );
          CPPUNIT_TEST (testIsZeroFarFromMembrane<stencil::TwoPoint> );
          CPPUNIT_TEST (testIsSymmetric<stencil::FourPoint> );
          CPPUNIT_TEST (testIsSymmetric<stencil::CosineApprox> );
          CPPUNIT_TEST (testIsSymmetric<stencil::ThreePoint> );
          CPPUNIT_TEST (testIsSymmetric<stencil::TwoPoint> );
          CPPUNIT_TEST (testIsIncreasing<stencil::FourPoint> );
          CPPUNIT_TEST (testIsIncreasing<stencil::CosineApprox> );
          CPPUNIT_TEST (testIsIncreasing<stencil::ThreePoint> );
          CPPUNIT_TEST (testIsIncreasing<stencil::TwoPoint> );
          CPPUNIT_TEST (testIsLinear<stencil::FourPoint> );
          // CPPUNIT_TEST (testIsLinear<stencil::CosineApprox>);
          CPPUNIT_TEST (testIsLinear<stencil::ThreePoint> );
          CPPUNIT_TEST (testIsLinear<stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> Kernel;

        public:
          void setUp();
          void setUpForces();

          // Moves fixture so barycenter is at given position and computes spread
          // force at origin
          template<class STENCIL>
          LatticeForceVector force_at_center(LatticePosition const &position);
          template<class STENCIL> void testIsZeroFarFromMembrane();
          template<class STENCIL> void testIsSymmetric();
          template<class STENCIL> void testIsIncreasing();
          template<class STENCIL> void testIsLinear();

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

      template<class STENCIL>
      LatticeForceVector CellForceSpreadTests::force_at_center(LatticePosition const &position)
      {
        mesh += position - mesh.GetBarycenter();
        helpers::ZeroOutForces(latDat);
        details::spreadForce2Grid<details::SpreadForces, STENCIL>(std::shared_ptr<CellBase>(&mesh,
                                                                                            [](CellBase*)
                                                                                            {}),
                                                                  details::SpreadForces(forces,
                                                                                        *latDat));
	// TODO #759: is this truncation OK?
        return latDat->GetSite(LatticeVector{center}).GetForce();
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

      template<class STENCIL> void CellForceSpreadTests::testIsZeroFarFromMembrane()
      {
        // Very far away from pancake samosa
        auto const border = Dimensionless(STENCIL::GetRange()) / 2e0;
        LatticeForceVector const faraway = force_at_center<STENCIL>(center
            + LatticePosition(0, 0, border + 1e0));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, faraway.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, faraway.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, faraway.z, 1e-8);

        auto const justTooFar = force_at_center<STENCIL>(center + LatticePosition(0, 0, border));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, justTooFar.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, justTooFar.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, justTooFar.z, 1e-8);

        auto const justInside = force_at_center<STENCIL>(center
            + LatticePosition(0, 0, border - 1e-4));
        CPPUNIT_ASSERT(not helpers::is_zero(justInside));
        CPPUNIT_ASSERT(justInside[0] > 1e-8 and justInside[0] < 1e-4);
      }

      template<class STENCIL> void CellForceSpreadTests::testIsSymmetric()
      {
        size_t const N(10);

        for (size_t i(0); i < N; ++i)
        {
          auto const range = Dimensionless(STENCIL::GetRange() * 0.5);
          auto const disp = range * (1e0 - Dimensionless(i) / Dimensionless(N));
          LatticePosition const displacement(0, 0, disp);
          LatticeForceVector const left = force_at_center<STENCIL>(center + displacement);
          LatticeForceVector const right = force_at_center<STENCIL>(center + displacement);
          CPPUNIT_ASSERT(helpers::is_zero(left - right));
        }
      }

      template<class STENCIL> void CellForceSpreadTests::testIsIncreasing()
      {
        size_t const N(10);
        LatticeForceVector last(0, 0, 0);

        for (size_t i(0); i < N; ++i)
        {
          auto const range = Dimensionless(STENCIL::GetRange() * 0.5);
          auto const d = range * (1e0 - Dimensionless(i + 1) / Dimensionless(N));
          LatticePosition const displacement(0, 0, d);
          LatticeForceVector const current = force_at_center<STENCIL>(center + displacement);
          CPPUNIT_ASSERT(current[0] > last[0]);
          CPPUNIT_ASSERT(current[1] > last[1]);
          CPPUNIT_ASSERT(current[2] > last[2]);
          last = current;
        }
      }

      template<class STENCIL> void CellForceSpreadTests::testIsLinear()
      {
        size_t const N(5);
        mesh = Cell(refine(MeshData { mesh.GetVertices(), mesh.GetFacets() }, 4));
        setUpForces();
        // x0, x1 should be further than 2 from the edges
        // Only linear if samosa appears as infinite plane
        // with sufficiently dense vertices
        LatticePosition const x0(center[0], center[1] - 0.5, center[2] - 0.1);
        LatticePosition const x1(center[0], center[1] + 0.5, center[2] - 0.1);
        LatticeForceVector const v0(force_at_center<STENCIL>(x0)), v1(force_at_center<STENCIL>(x1));

        LatticeForceVector const a( (v1 - v0) / (direction.Dot(x1) - direction.Dot(x0)));

        auto const reltol = std::vector<double> { 1e-3, 1e-4, 1e-5 }.at(STENCIL::GetRange() - 2);
        for (size_t i(0); i < N; ++i)
        {
          LatticePosition const x = (x1 - x0) * (Dimensionless(i + 1) / Dimensionless(N + 2)) + x0;
          LatticeForceVector const expected(a * (direction.Dot(x) - direction.Dot(x0)) + v0);
          auto const actual = force_at_center<STENCIL>(x);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.x, actual.x, std::max(expected.x * reltol, 1e-8));
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.y, actual.y, std::max(expected.y * reltol, 1e-8));
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.z, actual.z, std::max(expected.z * reltol, 1e-8));
        }
      }


      CPPUNIT_TEST_SUITE_REGISTRATION (CellForceSpreadTests);
    }
  }
}

#endif  // ONCE
