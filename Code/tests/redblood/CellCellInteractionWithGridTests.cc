// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/Cell.h"
#include "redblood/CellCell.h"
#include "redblood/stencil.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    namespace {
      LatticeDistance constexpr cutoff = 5.0;
      LatticeDistance constexpr halo = 2.0;

      auto const approx_zero = ApproxVector<double>{0}.Margin(1e-8);
    }

    // Work around the Catch quirk that TEMPLATE_TEST_CASE_METHOD
    // applies its type argument to the fixture, not only to the test
    // body (a la TEMPLATE_TEST_CASE).
    template<typename T>
    struct BigBox : public helpers::FourCubeBasedTestFixture<32> {
    };

    TEMPLATE_TEST_CASE_METHOD(BigBox,
			      "CellCellInteractionWithGridTests",
			      "[redblood]",
			      stencil::FourPoint, stencil::CosineApprox, stencil::ThreePoint, stencil::TwoPoint) {
      using STENCIL = TestType;

      auto cells = TwoPancakeSamosas<>(cutoff);

      // Place two nodes close enough for interactions
      LatticePosition const n0 { 15 - 0.1, 15.5, 15.5 };
      LatticePosition const n1 { 15 + 0.1, 15.5, 15.5 };
      (*cells.begin())->GetVertices().front() = n0;
      (*std::next(cells.begin()))->GetVertices().front() = n1;

      // Set forces to zero
      helpers::ZeroOutFOld(this->latDat);

      // Finds pairs, computes interaction, spread forces to lattice
      addCell2CellInteractions<STENCIL>(DivideConquerCells(cells, cutoff, halo),
					Node2NodeForce(1.0, halo),
					*this->latDat);

      // By symmetry, there are no forces on the lattice points equidistant from
      // the nodes
      REQUIRE(approx_zero == this->latDat->GetSite(15, 15, 15).GetForce());
      REQUIRE(approx_zero == this->latDat->GetSite(15, 14, 14).GetForce());
      REQUIRE(approx_zero == this->latDat->GetSite(15, 16, 16).GetForce());

      // There are non-zero opposite forces on the following nodes
      size_t delta(1);
      for (; 2 * delta <= STENCIL::GetRange(); ++delta) {
	REQUIRE(this->latDat->GetSite(14, 15, 15).GetForce().GetMagnitudeSquared() > 1e-8);
	auto const left = this->latDat->GetSite(15 + delta, 15, 15).GetForce();
	auto const right = this->latDat->GetSite(15 - delta, 15, 15).GetForce();
	REQUIRE(ApproxV(left) == -right);
	// The forces at (14, 15, 15) should be  in direction (-1, 0, 0)
	REQUIRE(right.Dot( util::Vector3D<double>{-1, 0, 0}) == Approx(std::abs(right.x)).margin(1e-8));
      }
      // This node is too far away
      REQUIRE(approx_zero == this->latDat->GetSite(15 + delta, 15, 15).GetForce());
      REQUIRE(approx_zero == this->latDat->GetSite(15 - delta, 15, 15).GetForce());
    }
  }
}
