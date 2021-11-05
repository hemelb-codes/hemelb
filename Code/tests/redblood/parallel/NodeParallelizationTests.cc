// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <algorithm>

#include <catch2/catch.hpp>

#include "tests/redblood/Fixtures.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    using namespace redblood::parallel;

    TEST_CASE("NodeParallelizationTests", "[redblood]") {
      using V = NodeCharacterizer::Process2NodesMap::mapped_type;
      auto approx = Approx(0).margin(1e-8);

      SECTION("testProperties") {
        NodeCharacterizer const nc( { { 0, V { 0, 2 } }, { 1, V { 1, 2 } }, { 2, V { } } });
        REQUIRE(nc.IsMidDomain(0) == true);
        REQUIRE(nc.IsMidDomain(1) == true);
        REQUIRE(nc.IsMidDomain(2) == false);
        REQUIRE(nc.IsBoundary(0) == false);
        REQUIRE(nc.IsBoundary(1) == false);
        REQUIRE(nc.IsBoundary(2) == true);
      }

      SECTION("testConstruction") {
        auto func = [](LatticePosition const &position) -> std::set<proc_t>
        {
          bool isZero = std::abs(position.x) < 1e1 and std::abs(position.y) < 1e1
          and std::abs(position.z) < 1e1;
          bool isOne = std::abs(position.x - 15e0) < 1e1 or not isZero;
          if(isZero and isOne)
          return
          { 0, 1};
          return
          { isZero? 0: 1};
        };

        NodeCharacterizer nc(func, { LatticePosition(0, 0, 0),
                                     LatticePosition(1e2, 0, 0),
                                     LatticePosition(7e0, 0, 0) });

        REQUIRE(nc[0].count(0));
        REQUIRE(nc[0].count(2));
        REQUIRE(nc[1].count(1));
        REQUIRE(nc[1].count(2));
      }

      SECTION("testReduceFrom") {
        std::vector<LatticePosition> reduced(3, { 0, 0, 0 });
        std::vector<LatticePosition> const incomming0 = { { 1, 0, 0 }, { 0, 2, 0 } };
        std::vector<LatticePosition> const incomming1 = { { 0, 1, 1 }, { 0, 0, 1 } };

        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        nc.ReduceFrom(reduced, 0, incomming0);
        std::vector<LatticePosition> const expected0 = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 0 } };
        for (auto const item : util::czip(expected0, reduced))
        {
          REQUIRE(std::get<0>(item).x == approx(std::get<1>(item).x));
          REQUIRE(std::get<0>(item).y == approx(std::get<1>(item).y));
          REQUIRE(std::get<0>(item).z == approx(std::get<1>(item).z));
        }

        nc.ReduceFrom(reduced, 1, incomming1);
        std::vector<LatticePosition> const expected1 = { { 1, 1, 1 }, { 0, 2, 0 }, { 0, 0, 1 } };
        for (auto const item : util::czip(expected1, reduced))
        {
          REQUIRE(std::get<0>(item).x == approx(std::get<1>(item).x));
          REQUIRE(std::get<0>(item).y == approx(std::get<1>(item).y));
          REQUIRE(std::get<0>(item).z == approx(std::get<1>(item).z));
        }

        nc.ReduceFrom(reduced, 2, { });
        for (auto const item : util::czip(expected1, reduced))
        {
          REQUIRE(std::get<0>(item).x == approx(std::get<1>(item).x));
          REQUIRE(std::get<0>(item).y == approx(std::get<1>(item).y));
          REQUIRE(std::get<0>(item).z == approx(std::get<1>(item).z));
        }
      }

      SECTION("testReduceFromAll") {
        std::vector<LatticePosition> reduced(3, { 0, 0, 0 });

        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        nc.ReduceFrom(reduced, { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 1, 1 }, { 0, 0, 1 } });
        std::vector<LatticePosition> const expected = { { 1, 1, 1 }, { 0, 2, 0 }, { 0, 0, 1 } };
        for (auto const item : util::czip(expected, reduced))
        {
          REQUIRE(std::get<0>(item).x == approx(std::get<1>(item).x));
          REQUIRE(std::get<0>(item).y == approx(std::get<1>(item).y));
          REQUIRE(std::get<0>(item).z == approx(std::get<1>(item).z));
        }
      }

      SECTION("testSpreadTo") {
        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        auto const actual = nc.SpreadTo( { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } });
        REQUIRE(size_t(3) == actual.first.size());
        REQUIRE(size_t(2) == actual.first[0]);
        REQUIRE(size_t(2) == actual.first[1]);
        REQUIRE(size_t(0) == actual.first[2]);
        std::vector<LatticePosition> const expected =
            { { 1, 0, 0 }, { 0, 2, 0 }, { 1, 0, 0 }, { 0, 0, 3 } };
        REQUIRE(expected.size() == actual.second.size());
        for (auto const item : util::czip(expected, actual.second))
        {
          REQUIRE(std::get<0>(item).x == approx(std::get<1>(item).x));
          REQUIRE(std::get<0>(item).y == approx(std::get<1>(item).y));
          REQUIRE(std::get<0>(item).z == approx(std::get<1>(item).z));
        }
      }

    }
  }
}
