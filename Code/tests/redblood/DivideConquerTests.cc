// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iterator>
#include <catch2/catch.hpp>

#include "redblood/DivideConquer.h"

namespace hemelb
{
  namespace tests
  {
    //namespace redblood
    TEST_CASE("DivideAndConquerTests", "[redblood]") {
      using DnC = redblood::DivideConquer<int>;

      SECTION("testDowngradeKey") {
        LatticeDistance const cutoff = 5e0;
        DnC dnc(cutoff);

        size_t const N = 5;
        LatticePosition const inputs[N] = { LatticePosition(2.5, 1.1, 3.3), LatticePosition(5.5,
                                                                                            1.1,
                                                                                            3.3),
                                            LatticePosition(5.5, -5.1, -3.3), LatticePosition(5.5,
                                                                                              -5.1,
                                                                                              10.3),
                                            LatticePosition(5.000000000001, -5.1, 10.3) };
        LatticeVector const expected[N] = { LatticeVector(0, 0, 0),
                                            LatticeVector(1, 0, 0),
                                            LatticeVector(1, -2, -1),
                                            LatticeVector(1, -2, 2),
                                            LatticeVector(1, -2, 2) };

        for (size_t i(0); i < N; ++i)
	  REQUIRE(LatticeVector::Zero() == dnc.DowngradeKey(inputs[i]) - expected[i]);
      }

      SECTION("Checks downgrading does not occur if type is already a key", "NoDowngradeKey") {
        LatticeDistance const cutoff = 5e0;
        DnC dnc(cutoff);

        LatticeVector const key(10, 5, 2);
        REQUIRE(dnc.DowngradeKey(key) == key);
      }

      SECTION("testAddToBox") {
        LatticeDistance const cutoff = 5e0;
        DnC dnc(cutoff);

        typedef DnC::iterator iterator;
        iterator const i_inserted = dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
        REQUIRE(dnc.size() == 1);
        REQUIRE(i_inserted->first == LatticeVector(-1, 0, 1));
        REQUIRE(i_inserted->second == 2);

        // Adds exact same item -> two separate copies since this is a multimap
        iterator const i_other = dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
        REQUIRE(dnc.size() == 2);
        REQUIRE(i_other->first == LatticeVector(-1, 0, 1));
        REQUIRE(i_other->second == 2);
      }

      SECTION("testAddToBoxAsKey") {
        LatticeDistance const cutoff = 5e0;
        DnC dnc(cutoff);

        typedef DnC::iterator iterator;
        LatticeVector const key(10, 5, 2);
        iterator const i_inserted = dnc.insert(key, 2);
        REQUIRE(dnc.size() == 1);
        REQUIRE(i_inserted->first == key);
        REQUIRE(i_inserted->second == 2);
      }

      SECTION("testBoxRange") {
        LatticeDistance const cutoff = 5e0;
        DnC dnc(cutoff);
        dnc.insert(LatticePosition(-3.5, 0.1, 5.1), 2);
        dnc.insert(LatticePosition(-3.6, 0.2, 6.1), 4);
        dnc.insert(LatticePosition(0, 0.1, 5.1), 2);
        dnc.insert(LatticePosition(0, 0.3, 5.1), 4);
        REQUIRE(dnc.size() == 4);

        // Checks we can access range directly using box indices
        DnC::const_range asInt = dnc.equal_range(LatticeVector(-1, 0, 1));
        REQUIRE(std::distance(asInt.first, asInt.second) == 2);
        REQUIRE(asInt.first != asInt.second);
        REQUIRE(asInt.first->first == LatticeVector(-1, 0, 1));
        // No order guarantee until C++11
        REQUIRE((asInt.first->second == 2 or asInt.first->second == 4));
        DnC::const_iterator i_other = asInt.first;
        ++i_other;
        REQUIRE(i_other->first == asInt.first->first);
        REQUIRE(i_other->second != asInt.first->second);
        REQUIRE((i_other->second == 2 or i_other->second == 4));

        // Checks we can access range using position
        DnC::const_range asFloat = dnc.equal_range(LatticePosition(-3.5, 0.1, 5.1));
        REQUIRE(asInt.first == asFloat.first);
        REQUIRE(asInt.second == asFloat.second);

        // Checks empty range
        DnC::const_range empty = dnc.equal_range(LatticeVector(-10, 0, 1));
        REQUIRE(empty.first == empty.second);
      }

      SECTION("testBoxRangeAsKey") {
        LatticeDistance const cutoff = 5e0;

        DnC dnc(cutoff);
        LatticePosition const position(5.2, -3.3, 0.1);
        LatticeVector const key(1, -1, 0);

        dnc.insert(key, 2);
        dnc.insert(position, 3);

        // Checks we can access range directly using box indices
        DnC::range const range = dnc.equal_range(key);
        DnC::const_range const crange = range;
        REQUIRE(range == dnc.equal_range(position));
        REQUIRE(crange == const_cast<DnC const&>(dnc).equal_range(position));
        REQUIRE(crange == const_cast<DnC const&>(dnc).equal_range(key));
      }

    }
  }
}
