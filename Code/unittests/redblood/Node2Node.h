//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_NODE2NODE_H
#define HEMELB_UNITTESTS_REDBLOOD_NODE2NODE_H

#include <cppunit/TestFixture.h>
#include "redblood/Node2Node.h"

namespace hemelb { namespace unittests {

class Node2NodeTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(Node2NodeTests);
      CPPUNIT_TEST(testNode2NodeForce);
    CPPUNIT_TEST_SUITE_END();

public:
    void testNode2NodeForce() {
      PhysicalDistance const cutoff = 2e0 * 2e0; // Squared
      CPPUNIT_ASSERT_DOUBLES_EQUAL(
          redblood::node2NodeForce(cutoff, 1e0, cutoff, 1),
          0e0, 1e-12
      );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(
          redblood::node2NodeForce(1.1 * cutoff, 1e0, cutoff, 1),
          0e0, 1e-12
      );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(
          redblood::node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1),
          -1.0 / cutoff * (1.0 / 0.9 - 1.e0), 1e-12
      );
      CPPUNIT_ASSERT_DOUBLES_EQUAL(
          redblood::node2NodeForce(0.9 * cutoff, 1.1, cutoff, 1),
          1.1 * redblood::node2NodeForce(0.9 * cutoff, 1e0, cutoff, 1), 1e-12
      );

      // Test direction
      LatticePosition const direction = LatticePosition(1, 2, 3).Normalise();
      LatticeForceVector const force = redblood::node2NodeForce(
          direction, 1, direction.GetMagnitudeSquared() * 1.1
      );
      CPPUNIT_ASSERT(helpers::is_zero(
          force + direction * force.GetMagnitude()
      ));

      LatticePosition const A(0, 0, 0);
      LatticePosition const B(std::sqrt(0.9*cutoff), 0, 0);
      CPPUNIT_ASSERT(helpers::is_zero(
            redblood::node2NodeForce(B - A, 1.0, cutoff)
            + B.GetNormalised() * (1.0 / cutoff * (1.0 / 0.9 - 1.e0))
      ));
      CPPUNIT_ASSERT(helpers::is_zero(
            redblood::node2NodeForce(A, B, 1.0, cutoff)
            + (B - A).GetNormalised() * (1.0 / cutoff * (1.0 / 0.9 - 1.e0))
      ));
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(Node2NodeTests);
}}

#endif // ONCE

