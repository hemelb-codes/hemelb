// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_NODEPARALLELIZATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_NODEPARALLELIZATIONTESTS_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "util/Iterator.h"
#include <algorithm>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
      class NodeParallelizationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (NodeParallelizationTests);
          CPPUNIT_TEST (testProperties);
          CPPUNIT_TEST (testConstruction);
          CPPUNIT_TEST (testReduceFrom);
          CPPUNIT_TEST (testReduceFromAll);
          CPPUNIT_TEST (testSpreadTo);CPPUNIT_TEST_SUITE_END();

        public:

          void testProperties();
          void testConstruction();
          void testReduceFrom();
          void testReduceFromAll();
          void testSpreadTo();
      };

      void NodeParallelizationTests::testProperties()
      {
        typedef NodeCharacterizer::Process2NodesMap::mapped_type V;
        NodeCharacterizer const nc( { { 0, V { 0, 2 } }, { 1, V { 1, 2 } }, { 2, V { } } });
        CPPUNIT_ASSERT_EQUAL(nc.IsMidDomain(0), true);
        CPPUNIT_ASSERT_EQUAL(nc.IsMidDomain(1), true);
        CPPUNIT_ASSERT_EQUAL(nc.IsMidDomain(2), false);
        CPPUNIT_ASSERT_EQUAL(nc.IsBoundary(0), false);
        CPPUNIT_ASSERT_EQUAL(nc.IsBoundary(1), false);
        CPPUNIT_ASSERT_EQUAL(nc.IsBoundary(2), true);
      }

      void NodeParallelizationTests::testConstruction()
      {
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

        CPPUNIT_ASSERT(nc[0].count(0));
        CPPUNIT_ASSERT(nc[0].count(2));
        CPPUNIT_ASSERT(nc[1].count(1));
        CPPUNIT_ASSERT(nc[1].count(2));
      }

      void NodeParallelizationTests::testReduceFrom()
      {
        std::vector<LatticePosition> reduced(3, { 0, 0, 0 });
        std::vector<LatticePosition> const incomming0 = { { 1, 0, 0 }, { 0, 2, 0 } };
        std::vector<LatticePosition> const incomming1 = { { 0, 1, 1 }, { 0, 0, 1 } };
        typedef NodeCharacterizer::Process2NodesMap::mapped_type V;
        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        nc.ReduceFrom(reduced, 0, incomming0);
        std::vector<LatticePosition> const expected0 = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 0 } };
        for (auto const item : util::czip(expected0, reduced))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }

        nc.ReduceFrom(reduced, 1, incomming1);
        std::vector<LatticePosition> const expected1 = { { 1, 1, 1 }, { 0, 2, 0 }, { 0, 0, 1 } };
        for (auto const item : util::czip(expected1, reduced))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }

        nc.ReduceFrom(reduced, 2, { });
        for (auto const item : util::czip(expected1, reduced))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }
      }

      void NodeParallelizationTests::testReduceFromAll()
      {
        std::vector<LatticePosition> reduced(3, { 0, 0, 0 });
        typedef NodeCharacterizer::Process2NodesMap::mapped_type V;
        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        nc.ReduceFrom(reduced, { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 1, 1 }, { 0, 0, 1 } });
        std::vector<LatticePosition> const expected = { { 1, 1, 1 }, { 0, 2, 0 }, { 0, 0, 1 } };
        for (auto const item : util::czip(expected, reduced))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }
      }

      void NodeParallelizationTests::testSpreadTo()
      {
        typedef NodeCharacterizer::Process2NodesMap::mapped_type V;
        NodeCharacterizer const nc( { { 0, V { 0, 1 } }, { 1, V { 0, 2 } }, { 2, V { } } });

        auto const actual = nc.SpreadTo( { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 3 } });
        CPPUNIT_ASSERT_EQUAL(size_t(3), actual.first.size());
        CPPUNIT_ASSERT_EQUAL(size_t(2), actual.first[0]);
        CPPUNIT_ASSERT_EQUAL(size_t(2), actual.first[1]);
        CPPUNIT_ASSERT_EQUAL(size_t(0), actual.first[2]);
        std::vector<LatticePosition> const expected =
            { { 1, 0, 0 }, { 0, 2, 0 }, { 1, 0, 0 }, { 0, 0, 3 } };
        CPPUNIT_ASSERT_EQUAL(expected.size(), actual.second.size());
        for (auto const item : util::czip(expected, actual.second))
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (NodeParallelizationTests);
    }
  }
}

#endif  // ONCE
