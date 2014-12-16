//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_DIVIDE_AND_CONQUER_H
#define HEMELB_UNITTESTS_REDBLOOD_DIVIDE_AND_CONQUER_H

#include <iterator>
#include <cppunit/TestFixture.h>
#include "redblood/DivideConquer.h"
#include "unittests/helpers/Comparisons.h"

namespace hemelb { namespace unittests { namespace redblood {

class DivideAndConquerTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(DivideAndConquerTests);
      CPPUNIT_TEST(testDowngradeKey);
      CPPUNIT_TEST(testNoDowngradeKey);
      CPPUNIT_TEST(testAddToBox);
      CPPUNIT_TEST(testAddToBoxAsKey);
      CPPUNIT_TEST(testBoxRange);
      CPPUNIT_TEST(testBoxRangeAsKey);
    CPPUNIT_TEST_SUITE_END();

    typedef DivideConquer<int> DnC;
public:
    void testDowngradeKey();
    void testAddToBox();
    void testBoxRange();
    // Checks downgrading does not occur if type is already a key
    void testNoDowngradeKey();
    void testAddToBoxAsKey();
    void testBoxRangeAsKey();
};


void DivideAndConquerTests :: testDowngradeKey() {
  PhysicalDistance const cutoff = 5e0;
  DnC dnc(cutoff);

  size_t const N=5;
  LatticePosition const inputs[N] = {
    LatticePosition(2.5, 1.1, 3.3),
    LatticePosition(5.5, 1.1, 3.3),
    LatticePosition(5.5, -5.1, -3.3),
    LatticePosition(5.5, -5.1, 10.3),
    LatticePosition(5.000000000001, -5.1, 10.3)
  };
  LatticeVector const expected[N] = {
    LatticeVector(0, 0, 0),
    LatticeVector(1, 0, 0),
    LatticeVector(1, -2, -1),
    LatticeVector(1, -2, 2),
    LatticeVector(1, -2, 2)
  };
  for(size_t i(0); i < N; ++i)
    CPPUNIT_ASSERT(helpers::is_zero(
      dnc.DowngradeKey(inputs[i]) - expected[i]
    ));
}

void DivideAndConquerTests :: testNoDowngradeKey() {
  PhysicalDistance const cutoff = 5e0;
  DnC dnc(cutoff);

  LatticeVector const key(10, 5, 2);
  CPPUNIT_ASSERT(dnc.DowngradeKey(key) == key);
}

void DivideAndConquerTests :: testAddToBox() {
  PhysicalDistance const cutoff = 5e0;
  DnC dnc(cutoff);

  typedef DnC::iterator iterator;
  iterator const i_inserted = dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
  CPPUNIT_ASSERT(dnc.size() == 1);
  CPPUNIT_ASSERT(i_inserted->first == LatticeVector(-1, 0, 1));
  CPPUNIT_ASSERT(i_inserted->second == 2);

  // Adds exact same item -> two separate copies since this is a multimap
  iterator const i_other = dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
  CPPUNIT_ASSERT(dnc.size() == 2);
  CPPUNIT_ASSERT(i_other->first == LatticeVector(-1, 0, 1));
  CPPUNIT_ASSERT(i_other->second == 2);
}

void DivideAndConquerTests :: testAddToBoxAsKey() {
  PhysicalDistance const cutoff = 5e0;
  DnC dnc(cutoff);

  typedef DnC::iterator iterator;
  LatticeVector const key(10, 5, 2);
  iterator const i_inserted = dnc.insert(key, 2);
  CPPUNIT_ASSERT(dnc.size() == 1);
  CPPUNIT_ASSERT(i_inserted->first == key);
  CPPUNIT_ASSERT(i_inserted->second == 2);
}

void DivideAndConquerTests :: testBoxRange() {
  PhysicalDistance const cutoff = 5e0;
  DnC dnc(cutoff);
  dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
  dnc.insert(LatticePosition(-3.6, 0.2, 6.1), 4);
  dnc.insert(LatticePosition(0, 0.1, 5.1), 2);
  dnc.insert(LatticePosition(0, 0.3, 5.1), 4);
  CPPUNIT_ASSERT(dnc.size() == 4);

  // Checks we can access range directly using box indices
  DnC::const_range asInt = dnc.equal_range(LatticeVector(-1, 0, 1));
  CPPUNIT_ASSERT(std::distance(asInt.first, asInt.second) == 2);
  CPPUNIT_ASSERT(asInt.first != asInt.second);
  CPPUNIT_ASSERT(asInt.first->first == LatticeVector(-1, 0, 1));
  // No order guarantee until C++11
  CPPUNIT_ASSERT(asInt.first->second == 2 or asInt.first->second == 4);
  DnC::const_iterator i_other = asInt.first; ++i_other;
  CPPUNIT_ASSERT(i_other->first == asInt.first->first);
  CPPUNIT_ASSERT(i_other->second != asInt.first->second);
  CPPUNIT_ASSERT(i_other->second == 2 or i_other->second == 4);


  // Checks we can access range using position
  DnC::const_range asFloat = dnc.equal_range(LatticePosition(-3.5, 0.1, 5.1));
  CPPUNIT_ASSERT(asInt.first == asFloat.first);
  CPPUNIT_ASSERT(asInt.second == asFloat.second);

  // Checks empty range
  DnC::const_range empty = dnc.equal_range(LatticeVector(-10, 0, 1));
  CPPUNIT_ASSERT(empty.first == empty.second);
}

void DivideAndConquerTests :: testBoxRangeAsKey() {
  PhysicalDistance const cutoff = 5e0;

  DnC dnc(cutoff);
  LatticePosition const position(5.2, -3.3, 0.1);
  LatticeVector const key(1, -1, 0);

  dnc.insert(key, 2);
  dnc.insert(position, 3);

  // Checks we can access range directly using box indices
  DnC::range const range = dnc.equal_range(key);
  DnC::const_range const crange = range;
  CPPUNIT_ASSERT(range == dnc.equal_range(position));
  CPPUNIT_ASSERT(crange == const_cast<DnC const&>(dnc).equal_range(position));
  CPPUNIT_ASSERT(crange == const_cast<DnC const&>(dnc).equal_range(key));
}


CPPUNIT_TEST_SUITE_REGISTRATION(DivideAndConquerTests);
}}}

#endif // ONCE

