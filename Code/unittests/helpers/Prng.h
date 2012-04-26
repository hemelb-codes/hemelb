#ifndef HEMELB_UNITTESTS_HELPERS_PRNG_H
#define HEMELB_UNITTESTS_HELPERS_PRNG_H

#include <limits>

namespace hemelb
{

  namespace unittests
  {

    namespace helpers
    {
      /**
       * Need a source of data; use a linear congruential PRNG as no
       * requirement for it to be "good". Implementation from Numerical
       * Recipes. This implementation must give identical output to the
       * Python implementation in MakeDummyExtraction.py (given the same seed)
       *
       */
      class Prng
      {
        private:
          /**
           * Multiplier
           */
          enum
          {
            a = 1664525
          //!< a
          };
          /**
           * Increment
           */
          enum
          {
            c = 1013904223
          //!< c
          };

          uint32_t state;

        public:
          Prng(uint32_t seed) :
            state(seed)
          {

          }
          /**
           * Get a random int by
           * State[n+1] = (a State[n] + c) mod m
           * Since we have m == 2**32 the overflow does the
           * modulus operation for us.
           */
          uint32_t rand()
          {
            state = (a * state + c);
            return state;
          }

          float uniform()
          {
            return rand() / float(std::numeric_limits<uint32_t>::max());
          }
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_PRNG_H
