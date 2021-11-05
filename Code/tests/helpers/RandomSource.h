// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_RANDOMSOURCE_H
#define HEMELB_TESTS_HELPERS_RANDOMSOURCE_H

#include <limits>
#include <cstdint>

namespace hemelb
{

  namespace tests
  {

    namespace helpers
    {
      // Need a source of data; use a linear congruential
      // pseudorandom number generator as there is no requirement for
      // it to be "good".
      //
      // Implementation from Numerical Recipes. This implementation must give
      // identical output to the Python implementation in
      // MakeDummyExtraction.py (given the same seed)
      class RandomSource
      {
      private:
	// Multiplier
	static constexpr std::uint32_t a = 1664525;
	// Increment
	static constexpr std::uint32_t c = 1013904223;

	std::uint32_t state;

      public:
	inline RandomSource(std::uint32_t seed) :  state(seed)
	{
	}
	// Get a random int by
	// State[n+1] = (a State[n] + c) mod m
	// Since we have m == 2**32 the overflow does the
	// modulus operation for us.
	inline std::uint32_t rand() {
	  state = (a * state + c);
	  return state;
	}

	inline float uniform() {
	  return rand() / float(std::numeric_limits<uint32_t>::max());
	}
      };

    }
  }
}

#endif // HEMELB_TESTS_HELPERS_RANDOMSOURCE_H
