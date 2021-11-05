// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "tests/helpers/RandomSource.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      TEST_CASE("RandomSource") {
	SECTION("Ensure that the RandomSource gives the same output as the Python implementation") {
	  // These are the output of the Python
	  const std::vector<unsigned> randoms = {
	    3274329173, 1951117744, 2806192463, 4020483682, 1912357465, 56543204, 2972185075,
	    3331854006, 1472870045, 2563458904, 4106658519, 3060130378, 822013217, 1947542540,
	    319459323, 2027504926, 2968480229, 4096102912, 3699596639, 3835005746, 90108137,
	    3707700532, 2475224835, 2699643014, 1883440685, 1051799976, 3734989031, 1289900314,
	    1002925489, 1905068892, 1327181579, 2913057006, 2358389621, 3684222544, 2684675439,
	    3085947906, 3595010937, 1933825668, 289239059, 2299541078, 4132305341, 861597688,
	    2195845879, 2731914218, 953307713, 4097574572, 888132123, 4277559486, 3423310085,
	    1581928096
	  };
	  REQUIRE(randoms.size() == 50);

	  RandomSource generator(1358);
	  for (auto& expected: randoms) {
	    auto rand = generator.rand();
	    REQUIRE(rand == expected);
	  }
	}
      }
    }
  }
}
