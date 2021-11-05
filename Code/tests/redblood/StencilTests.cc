// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "redblood/stencil.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    TEST_CASE("StencilTests", "[redblood]") {
      auto approx = Approx(0).margin(1e-8);

      SECTION("4point") {
	using redblood::stencil::fourPoint;
	int constexpr n = 11;
	Dimensionless const x[] = { 0,
				    1. - 1e-8,
				    1. + 1e-8,
				    -1. - 1e-8,
				    -1. + 1e-8,
				    2. - 1e-8,
				    2. + 1e-8,
				    -2. - 1e-8,
				    -2. + 1e-8,
				    2.1,
				    -2.1};
	static_assert(sizeof(x) == sizeof(Dimensionless)*n, "Array size mismatch");
	Dimensionless const expected[] = { 0.5, 0.25, 0.25, 0.25, 0.25, 0., 0., 0., 0., 0., 0. };
	static_assert(sizeof(expected) == sizeof(Dimensionless)*n, "Array size mismatch");

	for (size_t i = 0; i < n; ++i) {
	  Dimensionless const actual = fourPoint(x[i]);
	  REQUIRE(actual == approx(expected[i]));
	}
      }

      SECTION("CosineApprox") {
	using redblood::stencil::cosineApprox;
	int constexpr n = 7;
	Dimensionless const x[] = { 0,
				    1. - 1e-8,
				    -1. + 1e-8,
				    2. - 1e-8,
				    -2. + 1e-8,
				    2.1,
				    -2.1
	};
	static_assert(sizeof(x) == sizeof(Dimensionless)*n, "Array size mismatch");
	Dimensionless const expected[] = { 0.5, 0.25, 0.25, 0., 0., 0., 0. };
	static_assert(sizeof(expected) == sizeof(Dimensionless)*n, "Array size mismatch");

	for (size_t i = 0; i < n; ++i) {
	  Dimensionless const actual = cosineApprox(x[i]);
	  REQUIRE(actual == approx(expected[i]));
	}
      }

      SECTION("3point") {
	using redblood::stencil::threePoint;
	int constexpr n = 11;
	Dimensionless const x[] = { 0, 0.5 - 1e-8, 0.5 + 1e-8, -0.5 - 1e-8, -0.5 + 1e-8, 1.5
				    - 1e-8,
				    1.5 + 1e-8, -1.5 - 1e-8, -1.5 + 1e-8, 1.6, -1.6
	};
	static_assert(sizeof(x) == sizeof(Dimensionless)*n, "Array size mismatch");
	Dimensionless expected[] = { 2. / 3., 0.5, 0.5, 0.5, 0.5, 0., 0., 0., 0., 0., 0. };
	static_assert(sizeof(expected) == sizeof(Dimensionless)*n, "Array size mismatch");

	for (size_t i = 0; i < n; ++i) {
	  Dimensionless const actual = threePoint(x[i]);
	  INFO(i);
	  REQUIRE(actual == approx(expected[i]));
	}
      }

      SECTION("2point") {
	using redblood::stencil::twoPoint;
	int constexpr n = 8;
	Dimensionless const x[] = { 0,
				    0.5,
				    1. - 1e-9,
				    1. + 1e-9,
				    -1. - 1e-9,
				    -1. + 1e-9,
				    1.1,
				    -1.1,
	};
	static_assert(sizeof(x) == sizeof(Dimensionless)*n, "Array size mismatch");
	Dimensionless expected[] = { 1., 0.5, 0., 0., 0., 0., 0., 0. };
	static_assert(sizeof(expected) == sizeof(Dimensionless)*n, "Array size mismatch");

	for (size_t i = 0; i < n; ++i) {
	  Dimensionless const actual = twoPoint(x[i]);
	  REQUIRE(actual == approx(expected[i]));
	}
      }
    }
  }
}
