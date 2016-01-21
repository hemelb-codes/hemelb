
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_RANDOMSOURCETESTS_H
#define HEMELB_UNITTESTS_HELPERS_RANDOMSOURCETESTS_H

#include <cppunit/TestFixture.h>

#include "unittests/helpers/RandomSource.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class RandomSourceTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE( RandomSourceTests);
          CPPUNIT_TEST( TestRand);CPPUNIT_TEST_SUITE_END();
        public:

          /**
           * Ensure that the RandomSource gives the same output as the Python implementation
           */
          void TestRand()
          {
            size_t nRandoms = 50;
            unsigned seed = 1358;
            // These are the output of the Python
            unsigned randoms[] = { 3274329173, 1951117744, 2806192463, 4020483682, 1912357465, 56543204, 2972185075,
                                   3331854006, 1472870045, 2563458904, 4106658519, 3060130378, 822013217, 1947542540,
                                   319459323, 2027504926, 2968480229, 4096102912, 3699596639, 3835005746, 90108137,
                                   3707700532, 2475224835, 2699643014, 1883440685, 1051799976, 3734989031, 1289900314,
                                   1002925489, 1905068892, 1327181579, 2913057006, 2358389621, 3684222544, 2684675439,
                                   3085947906, 3595010937, 1933825668, 289239059, 2299541078, 4132305341, 861597688,
                                   2195845879, 2731914218, 953307713, 4097574572, 888132123, 4277559486, 3423310085,
                                   1581928096 };

            RandomSource generator(seed);
            for (unsigned i = 0; i < nRandoms; ++i)
            {
              unsigned rand = generator.rand();
              CPPUNIT_ASSERT_EQUAL(randoms[i], rand);
            }
          }

      };

      CPPUNIT_TEST_SUITE_REGISTRATION( RandomSourceTests);
    }
  }
}
#endif // HEMELB_UNITTESTS_HELPERS_RANDOMSOURCETESTS_H
