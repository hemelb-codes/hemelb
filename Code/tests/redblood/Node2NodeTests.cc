// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/Node2Node.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace tests
  {
    TEST_CASE("Node2NodeTests", "[redblood]") {
      using namespace redblood;

      LatticeDistance const cutoff = 2e0;
      auto zero = Approx(0.0).margin(1e-12);
      REQUIRE(node2NodeForce(cutoff, 1e0, cutoff, 1) == zero);
      REQUIRE(node2NodeForce(1.1 * cutoff, 1e0, cutoff, 1) == zero);
      REQUIRE(node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1) == zero(-1.0 / cutoff * (1.0 / 0.9 - 1.e0)));
      REQUIRE(node2NodeForce(0.9 * cutoff, 1.1, cutoff, 1) == zero(1.1 * node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1)));

      // Test direction
      LatticePosition const direction = LatticePosition(1, 2, 3).Normalise();
      LatticeForceVector const force = node2NodeForce(direction,
						      1,
						      direction.GetMagnitude() * 1.1);
      REQUIRE(force == ApproxV(direction * (-force.GetMagnitude())));

      LatticePosition const A(0, 0, 0);
      LatticePosition const B(std::sqrt(0.9) * cutoff, 0, 0);
      LatticeForceVector const expected = -B.GetNormalised()
	* (1.0 / cutoff / cutoff * (1.0 / 0.9 - 1.e0));
      REQUIRE(node2NodeForce( (B - A).GetMagnitude(), 1.0, cutoff) == Approx(-expected.GetMagnitude()));
      REQUIRE(node2NodeForce(B - A, 1.0, cutoff) == ApproxV(expected));
      REQUIRE(node2NodeForce(A, B, 1.0, cutoff) == ApproxV(expected));
    }
  }
}
