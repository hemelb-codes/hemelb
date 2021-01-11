// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_NODE2NODETESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_NODE2NODETESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/Node2Node.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class Node2NodeTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (Node2NodeTests);
          CPPUNIT_TEST (testNode2NodeForce);CPPUNIT_TEST_SUITE_END();

        public:
          void testNode2NodeForce()
          {
            LatticeDistance const cutoff = 2e0;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(node2NodeForce(cutoff, 1e0, cutoff, 1), 0e0, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(node2NodeForce(1.1 * cutoff, 1e0, cutoff, 1), 0e0, 1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1),
                                         -1.0 / cutoff * (1.0 / 0.9 - 1.e0),
                                         1e-12);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(node2NodeForce(0.9 * cutoff, 1.1, cutoff, 1),
                                         1.1 * node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1),
                                         1e-12);

            // Test direction
            LatticePosition const direction = LatticePosition(1, 2, 3).Normalise();
            LatticeForceVector const force = node2NodeForce(direction,
                                                            1,
                                                            direction.GetMagnitude() * 1.1);
            CPPUNIT_ASSERT(helpers::is_zero(force + direction * force.GetMagnitude()));

            LatticePosition const A(0, 0, 0);
            LatticePosition const B(std::sqrt(0.9) * cutoff, 0, 0);
            LatticeForceVector const expected = -B.GetNormalised()
                * (1.0 / cutoff / cutoff * (1.0 / 0.9 - 1.e0));
            CPPUNIT_ASSERT(helpers::is_zero(node2NodeForce( (B - A).GetMagnitude(), 1.0, cutoff)
                + expected.GetMagnitude()));
            CPPUNIT_ASSERT(helpers::is_zero(node2NodeForce(B - A, 1.0, cutoff) - expected));
            CPPUNIT_ASSERT(helpers::is_zero(node2NodeForce(A, B, 1.0, cutoff) - expected));
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (Node2NodeTests);
    }
  }
}

#endif  // ONCE
