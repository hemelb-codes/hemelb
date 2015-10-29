//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_NODEPARALLELIZATION_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_NODEPARALLELIZATION_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/parallel/NodeCharacterizer.h"
// to unittest helper functions in in anonymous namespace
// #include "redblood/parallel/NodeCharacterizer.cc"
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
          CPPUNIT_TEST (testApplyMidDomain);
          CPPUNIT_TEST_SUITE_END();

        public:

          void testProperties();
          void testConstruction();
          void testApplyMidDomain();
      };

      void NodeParallelizationTests::testProperties()
      {
        NodeCharacterizer mm({{0}, {1}, {0, 1}});
        CPPUNIT_ASSERT_EQUAL(mm.IsMidDomain(0), true);
        CPPUNIT_ASSERT_EQUAL(mm.IsMidDomain(1), true);
        CPPUNIT_ASSERT_EQUAL(mm.IsMidDomain(2), false);
        CPPUNIT_ASSERT_EQUAL(mm.IsBoundary(0), false);
        CPPUNIT_ASSERT_EQUAL(mm.IsBoundary(1), false);
        CPPUNIT_ASSERT_EQUAL(mm.IsBoundary(2), true);
      }

      void NodeParallelizationTests::testConstruction()
      {
        auto func = [](LatticePosition const &position) -> NodeCharacterizer::ProcessorSet
        {
          bool isZero = std::abs(position.x) < 1e1 and std::abs(position.y) < 1e1
            and std::abs(position.z) < 1e1;
          bool isOne = std::abs(position.x - 15e0) < 1e1 or not isZero;
          if(isZero and isOne)
            return {0, 1};
          return {isZero? 0: 1};
        };

        NodeCharacterizer mm(func,
            {LatticePosition(0, 0, 0), LatticePosition(1e2, 0, 0), LatticePosition(7e0, 0, 0)});

        CPPUNIT_ASSERT(NodeCharacterizer::ProcessorSet{0} == mm[0]);
        CPPUNIT_ASSERT(NodeCharacterizer::ProcessorSet{1} == mm[1]);
        CPPUNIT_ASSERT(NodeCharacterizer::ProcessorSet({0, 1}) == mm[2]);
      }

      void NodeParallelizationTests::testApplyMidDomain()
      {
        auto const cell = std::make_shared<Cell>(pancakeSamosa());
        auto positionOnProc = [](LatticePosition const &position) -> std::set<proc_t>
        {
          bool isZero = position.GetMagnitude() < 1e-8;
          bool isSecond = (position - LatticePosition{1, 0, 1}).GetMagnitude() < 1e-8;
          return isZero ?
            std::set<proc_t>{0}: (isSecond ? std::set<proc_t>{1}: std::set<proc_t>{0, 1, 2});
        };
        MeshOwner const mm(positionOnProc, cell);
        size_t nbcalls(0);
        LatticePosition position(10, 10, 10);
        auto apply = [&position, &nbcalls](LatticePosition &input)
        {
          ++nbcalls;
          position = input;
        };

        // Finally, makes a call: there should be one node on proc 0
        mm.ApplyToMidDomainNodes(apply, 0);
        CPPUNIT_ASSERT_EQUAL(1ul, nbcalls);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, position.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, position.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, position.z, 1e-8);

        // there should be one node on proc 1
        mm.ApplyToMidDomainNodes(apply, 1);
        CPPUNIT_ASSERT_EQUAL(2ul, nbcalls);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1, position.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0, position.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1, position.z, 1e-8);

        // The node on proc 2 is shared (not mid-domain)
        mm.ApplyToMidDomainNodes(apply, 2);
        CPPUNIT_ASSERT_EQUAL(2ul, nbcalls);
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (NodeParallelizationTests);
    }
  }
}

#endif  // ONCE
