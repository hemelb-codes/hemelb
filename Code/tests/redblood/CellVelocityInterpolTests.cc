// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/kernels/LBGK.h"
#include "redblood/Facet.h"
#include "redblood/GridAndCell.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    using D3Q15 = lb::lattices::D3Q15;
    using Kernel = lb::kernels::LBGK<D3Q15>;

    class CellVelocityInterpolTests : public SquareDuctTetrahedronFixture
    {
    public:
      LatticePosition gradient;
      Dimensionless non_neg_pop;
      std::function<Dimensionless(LatticeVelocity const &)> linear, linear_inv;

      CellVelocityInterpolTests() : SquareDuctTetrahedronFixture{pancakeSamosa(0), 3}
      {
      }

      // Sets up a linear velocity profile
      void setupGradient(PhysicalVelocity const& GRADIENT) {
	std::tie(non_neg_pop, gradient, linear, linear_inv) =
	  helpers::makeLinearProfile(cubeSizeWithHalo, latDat, GRADIENT);
      }
      
    };

    // Checks fixture function do what they should do
    TEST_CASE_METHOD(CellVelocityInterpolTests, "testDistributionFixture", "[redblood]") {
      helpers::ZeroOutFOld(latDat);
      setupGradient(LatticeVelocity(2., 4., 6.));
      // Test assumes static pop at index == 0 as assumed by macro
      auto approx = [](double x) {
	return Approx(x).margin(1e-8);
      };
      auto zero = approx(0.0);
      REQUIRE(D3Q15::CX[0] == zero);
      REQUIRE(D3Q15::CY[0] == zero);
      REQUIRE(D3Q15::CZ[0] == zero);

      size_t constexpr N = 4;
      LatticeVector const a[N] = { LatticeVector(2, 4, 3),
				   LatticeVector(10, 16, 5),
				   LatticeVector(20, 3, 10),
				   LatticeVector(22, 8, 15) };

      for (size_t i = 0; i < N; ++i) {
	size_t const index = latDat->GetContiguousSiteId(a[i]);
	LatticePosition const pos(a[i][0], a[i][1], a[i][2]);
	REQUIRE(linear(pos) == approx(latDat->GetSite(index).GetFOld<D3Q15>()[0]));
	REQUIRE(linear_inv(pos) == approx(latDat->GetSite(index).GetFOld<D3Q15>()[1]));

	for (size_t j = 2; j < D3Q15::NUMVECTORS; ++j) {
	  REQUIRE(zero == latDat->GetSite(index).GetFOld<D3Q15>()[j]);
	}

	REQUIRE(approx(3.0 * non_neg_pop) ==
		latDat->GetSite(index).GetFOld<D3Q15>()[0] + latDat->GetSite(index).GetFOld<D3Q15>()[1]);
      }
    }

    // Becase Catch doesn't support template test cases + fixture,
    // only test cases over a templated fixture.
    template <typename T>
    class CellVelocityInterpolStencil : public CellVelocityInterpolTests {
    };

    using StencilTypes = std::tuple<stencil::FourPoint, stencil::ThreePoint, stencil::TwoPoint>;

    // Linear velociy profile on grid, perpendicular to the samosa
    // The interpolated velocities should be constant across the samosa
    TEMPLATE_LIST_TEST_CASE_METHOD(CellVelocityInterpolStencil,
				   "testLinearVelocityPerpendicularToPancakeSamosa",
				   "[redblood]",
				   StencilTypes) {
      using STENCIL = TestType;
      // direction perpendicular to plane
      helpers::ZeroOutFOld(this->latDat);
      this->setupGradient(Facet(*this->mesh.GetTemplateMesh().GetData(), 0).normal());

      // Perform interpolation
      std::vector<LatticePosition> displacements;
      // shared pointer with a fake deleter
      std::shared_ptr<CellBase> ptr_mesh(&this->mesh, [](CellBase*)
					 {});
      velocitiesOnMesh<Kernel, STENCIL>(ptr_mesh, *this->latDat, displacements);

      // Compute expected velocities
      using const_iterator = std::vector<LatticePosition>::const_iterator;
      const_iterator i_disp = displacements.begin();
      const_iterator const i_end = displacements.end();
      LatticePosition const expected(*i_disp);

      for (++i_disp; i_disp != i_end; ++i_disp) {
	REQUIRE(*i_disp == ApproxV(expected));
      }
    }

    // Linear velociy profile on grid, in-plane with the samosa
    // The interpolated velocities should evolve linearly with respect to the
    // input gradient.
    TEMPLATE_LIST_TEST_CASE_METHOD(CellVelocityInterpolStencil,
				   "testLinearVelocityInSamosaPlane",
				   "[redblood]",
				   StencilTypes) {
      using STENCIL = TestType;
      auto& mesh = this->mesh;
      // Figures out an in-plane direction
      helpers::ZeroOutFOld(this->latDat);
      Facet const shapeFacet(*mesh.GetTemplateMesh().GetData(), 0);
      LatticePosition const inplane(shapeFacet.edge(0) + shapeFacet.edge(1) * 0.5);
      this->setupGradient(inplane);

      // Perform interpolation
      std::vector<LatticePosition> displacements;
      std::shared_ptr<CellBase> ptr_mesh(&mesh, [](CellBase*)
					 {});
      velocitiesOnMesh<Kernel, STENCIL>(ptr_mesh, *this->latDat, displacements);

      // Computes what the interpolation should be
      using const_iterator = std::vector<LatticePosition>::const_iterator;
      auto const x0 = this->gradient.Dot(mesh.GetVertices()[0]);
      auto const x1 = this->gradient.Dot(mesh.GetVertices()[1]);
      LatticeVelocity const v0 = displacements[0], v1 = displacements[1];
      MeshData::Vertices::const_iterator i_vertex(mesh.GetVertices().begin() + 2);
      const_iterator i_disp = displacements.begin() + 2;
      const_iterator const i_end = displacements.end();

      for (; i_disp != i_end; ++i_disp, ++i_vertex) {
	LatticeVelocity const expected( (v0 - v1) * ( (i_vertex->Dot(this->gradient) - x1) / (x0 - x1))
					+ v1);

	REQUIRE(*i_disp == ApproxV(expected));
      }
    }

  }
}
