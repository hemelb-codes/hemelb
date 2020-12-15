// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLVELOCITYINTERPOLTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLVELOCITYINTERPOLTESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/GridAndCell.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      namespace
      {
        // Creates a mesh that is a single planar triangle
        // It still has two facets so that we can run fake forces on it.
        Mesh triangleMesh()
        {
          // Rotate samosa so it is in xy plane
          Mesh result = pancakeSamosa(0);
          result.GetData()->vertices[1] = LatticePosition(0, 1, 0);
          result.GetData()->vertices[2] = LatticePosition(std::sqrt(3.0) / 2.0, 0.5, 0);
          return result;
        }
      }

      class CellVelocityInterpolTests : public SquareDuctTetrahedronFixture
      {
          CPPUNIT_TEST_SUITE (CellVelocityInterpolTests);
          CPPUNIT_TEST (testDistributionFixture);
          CPPUNIT_TEST (testLinearVelocityPerpendicularToPancakeSamosa<stencil::FourPoint> );
          CPPUNIT_TEST (testLinearVelocityPerpendicularToPancakeSamosa<stencil::ThreePoint> );
          CPPUNIT_TEST (testLinearVelocityPerpendicularToPancakeSamosa<stencil::TwoPoint> );
          CPPUNIT_TEST (testLinearVelocityInSamosaPlane<stencil::FourPoint> );
          CPPUNIT_TEST (testLinearVelocityInSamosaPlane<stencil::ThreePoint> );
          CPPUNIT_TEST (testLinearVelocityInSamosaPlane<stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> Kernel;

        public:
          void setUp()
          {
            SquareDuctTetrahedronFixture::setUp();
          }

          // Checks fixture function do what they should do
          void testDistributionFixture();

          // Linear velociy profile on grid, perpendicular to the samosa
          // The interpolated velocities should be constant across the samosa
          template<class STENCIL> void testLinearVelocityPerpendicularToPancakeSamosa();
          // Linear velociy profile on grid, in-plane with the samosa
          // The interpolated velocities should evolve linearly with respect to the
          // input gradient.
          template<class STENCIL> void testLinearVelocityInSamosaPlane();

        protected:
          // Creates a mesh that is a single planar triangle
          // It still has two facets so that we can run fake forces on it.
          virtual Mesh initial_mesh() const
          {
            return pancakeSamosa(0);
          }
          virtual size_t refinement() const
          {
            return 3;
          }
      };

// Sets up a linear velocity profile
// Uses a macro so we can define a number of variables in one go
#define HEMELB_LINEAR_VELOCITY_PROFILE(GRADIENT)                             \
  LatticePosition gradient;                                                  \
  Dimensionless non_neg_pop;                                                 \
  std::function<Dimensionless(LatticeVelocity const &)> linear, linear_inv;  \
  std::tie(non_neg_pop, gradient, linear, linear_inv) =                      \
    helpers::makeLinearProfile(CubeSize(), latDat, GRADIENT);

      void CellVelocityInterpolTests::testDistributionFixture()
      {
        helpers::ZeroOutFOld(latDat);

        HEMELB_LINEAR_VELOCITY_PROFILE(LatticeVelocity(2., 4., 6.));
        // Test assumes static pop at index == 0 as assumed by macro
        CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CX[0], 0e0, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CY[0], 0e0, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(D3Q15::CZ[0], 0e0, 1e-8);

        size_t const N = 4;
        LatticeVector const a[N] = { LatticeVector(2, 4, 3),
                                     LatticeVector(10, 16, 5),
                                     LatticeVector(20, 3, 10),
                                     LatticeVector(22, 8, 15) };

        for (size_t i(0); i < N; ++i)
        {
          size_t const index = latDat->GetContiguousSiteId(a[i]);
          LatticePosition const pos(a[i][0], a[i][1], a[i][2]);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(linear(pos),
                                       latDat->GetSite(index).GetFOld<D3Q15>()[0],
                                       1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(linear_inv(pos),
                                       latDat->GetSite(index).GetFOld<D3Q15>()[1],
                                       1e-8);

          for (size_t j(2); j < D3Q15::NUMVECTORS; ++j)
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, latDat->GetSite(index).GetFOld<D3Q15>()[j], 1e-8);
          }

          CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0 * non_neg_pop,
                                       latDat->GetSite(index).GetFOld<D3Q15>()[0]
                                           + latDat->GetSite(index).GetFOld<D3Q15>()[1],
                                       1e-8);
        }
      }

      template<class STENCIL>
      void CellVelocityInterpolTests::testLinearVelocityPerpendicularToPancakeSamosa()
      {
        // direction perpendicular to plane
        helpers::ZeroOutFOld(latDat);
        HEMELB_LINEAR_VELOCITY_PROFILE(Facet(*mesh.GetTemplateMesh().GetData(), 0).normal());

        // Perform interpolation
        std::vector<LatticePosition> displacements;
        // shared pointer with a fake deleter
        std::shared_ptr<CellBase> ptr_mesh(&mesh, [](CellBase*)
        {});
        velocitiesOnMesh<Kernel, STENCIL>(ptr_mesh, *latDat, displacements);

        // Compute expected velocities
        typedef std::vector<LatticePosition>::const_iterator const_iterator;
        const_iterator i_disp = displacements.begin();
        const_iterator const i_end = displacements.end();
        LatticePosition const expected(*i_disp);

        for (++i_disp; i_disp != i_end; ++i_disp)
        {
          CPPUNIT_ASSERT(helpers::is_zero(*i_disp - expected));
        }
      }

      template<class STENCIL> void CellVelocityInterpolTests::testLinearVelocityInSamosaPlane()
      {
        // Figures out an in-plane direction
        helpers::ZeroOutFOld(latDat);
        Facet const shapeFacet(*mesh.GetTemplateMesh().GetData(), 0);
        LatticePosition const inplane(shapeFacet.edge(0) + shapeFacet.edge(1) * 0.5);
        HEMELB_LINEAR_VELOCITY_PROFILE(inplane);

        // Perform interpolation
        std::vector<LatticePosition> displacements;
        std::shared_ptr<CellBase> ptr_mesh(&mesh, [](CellBase*)
        {});
        velocitiesOnMesh<Kernel, STENCIL>(ptr_mesh, *latDat, displacements);

        // Computes what the interpolation should be
        typedef std::vector<LatticePosition>::const_iterator const_iterator;
        LatticeDistance const x0 = gradient.Dot(mesh.GetVertices()[0]), x1 =
            gradient.Dot(mesh.GetVertices()[1]);
        LatticeVelocity const v0 = displacements[0], v1 = displacements[1];
        MeshData::Vertices::const_iterator i_vertex(mesh.GetVertices().begin() + 2);
        const_iterator i_disp = displacements.begin() + 2;
        const_iterator const i_end = displacements.end();

        for (; i_disp != i_end; ++i_disp, ++i_vertex)
        {
          LatticeVelocity const expected( (v0 - v1) * ( (i_vertex->Dot(gradient) - x1) / (x0 - x1))
              + v1);
          CPPUNIT_ASSERT(helpers::is_zero(*i_disp - expected));
        }
      }
#undef HEMELB_LINEAR_VELOCITY_PROFILE

      CPPUNIT_TEST_SUITE_REGISTRATION (CellVelocityInterpolTests);
    }
  }
}

#endif  // ONCE
