// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_STENCILTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_STENCILTESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/stencil.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class StencilTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (StencilTests);
          CPPUNIT_TEST (test4point);
          CPPUNIT_TEST (test3point);
          CPPUNIT_TEST (test2point);
          CPPUNIT_TEST (testCosineApprox);CPPUNIT_TEST_SUITE_END();

        public:
          void test4point()
          {
            using stencil::fourPoint;
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
                                        -2.1,
                                        1e8 };
            Dimensionless expected[] = { 0.5, 0.25, 0.25, 0.25, 0.25, 0., 0., 0., 0., 0., 0. };

            for (size_t i(0); x[i] < 1e7; ++i)
            {
              Dimensionless const actual = fourPoint(x[i]);
              CPPUNIT_ASSERT(std::abs(actual - expected[i]) < 1e-8);
            }
          }
          void testCosineApprox()
          {
            using stencil::cosineApprox;
            Dimensionless const x[] = { 0,
                                        1. - 1e-8,
                                        -1. + 1e-8,
                                        2. - 1e-8,
                                        -2. + 1e-8,
                                        2.1,
                                        -2.1,
                                        1e8 // break value
                };
            Dimensionless expected[] = { 0.5, 0.25, 0.25, 0., 0., 0., 0. };

            for (size_t i(0); x[i] < 1e7; ++i)
            {
              Dimensionless const actual = cosineApprox(x[i]);
              CPPUNIT_ASSERT(std::abs(actual - expected[i]) < 1e-8);
            }
          }
          void test3point()
          {
            using stencil::threePoint;
            Dimensionless const x[] = { 0, 0.5 - 1e-8, 0.5 + 1e-8, -0.5 - 1e-8, -0.5 + 1e-8, 1.5
                                            - 1e-8,
                                        1.5 + 1e-8, -1.5 - 1e-8, -1.5 + 1e-8, 1.6, -1.6, 1e8 // break value
                };
            Dimensionless expected[] = { 2. / 3., 0.5, 0.5, 0.5, 0.5, 0., 0., 0., 0., 0., 0. };

            for (size_t i(0); x[i] < 1e7; ++i)
            {
              Dimensionless const actual = threePoint(x[i]);
              CPPUNIT_ASSERT(std::abs(actual - expected[i]) < 1e-8);
            }
          }
          void test2point()
          {
            using stencil::twoPoint;
            Dimensionless const x[] = { 0,
                                        0.5,
                                        1. - 1e-9,
                                        1. + 1e-9,
                                        -1. - 1e-9,
                                        -1. + 1e-9,
                                        1.1,
                                        -1.1,
                                        1e8 // break value
                };
            Dimensionless expected[] = { 1., 0.5, 0., 0., 0., 0., 0., 0. };

            for (size_t i(0); x[i] < 1e7; ++i)
            {
              Dimensionless const actual = twoPoint(x[i]);
              CPPUNIT_ASSERT(std::abs(actual - expected[i]) < 1e-8);
            }
          }
      };

      // empty line courtosy of cppunit
      CPPUNIT_TEST_SUITE_REGISTRATION (StencilTests);
    }
  }
}

#endif  // ONCE
