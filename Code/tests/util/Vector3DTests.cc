// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "util/Vector3D.h"


namespace hemelb
{
  namespace tests
  {

    TEST_CASE("Vector3D CastsInVector3DProduct") {
	const double dblMax = std::numeric_limits<double>::max();
	const unsigned uintMax = std::numeric_limits<unsigned>::max();

	SECTION("(int, float) -> float") {
	  util::Vector3D<int> foo(-1, 0, 1);
	  float bar = 0.8;
	  util::Vector3D<float> baz = foo * bar;

	  REQUIRE(-0.8f == baz[0]);
	  REQUIRE(0.0f == baz[1]);
	  REQUIRE(0.8f == baz[2]);
	}

	SECTION("(float, double) -> double"){
	  util::Vector3D<float> foo(0.0f, 1.0f, -1.0f);
	  double bar = dblMax;
	  util::Vector3D<double> baz = foo * bar;

	  REQUIRE(0.0 == baz[0]);
	  REQUIRE(dblMax == baz[1]);
	  REQUIRE(-dblMax == baz[2]);
	}

	SECTION("(int, unsigned) -> unsigned") {
	  util::Vector3D<int> foo(0, 2, 2);
	  unsigned bar = uintMax / 2u;
	  util::Vector3D<unsigned> baz = foo * bar;

	 REQUIRE(0u == baz[0]);
	 REQUIRE(uintMax == baz[1] + uintMax % 2);
	 REQUIRE(uintMax == baz[2] + uintMax % 2);
	}
    }

    TEST_CASE("Vector3D works at compile time") {
      using V = util::Vector3D<int>;

      constexpr auto z = V::Zero();
      STATIC_REQUIRE(z.x == 0);
      STATIC_REQUIRE(z.y == 0);
      STATIC_REQUIRE(z.z == 0);

      constexpr auto one = V::Ones();
      STATIC_REQUIRE(one.x == 1);
      STATIC_REQUIRE(one.y == 1);
      STATIC_REQUIRE(one.z == 1);

    }
  }
}
