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
      class MPINodeParallelizationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (MPINodeParallelizationTests);
          CPPUNIT_TEST(testNeighborhoodGraphCreation);
          CPPUNIT_TEST_SUITE_END();

        public:
          void testNeighborhoodGraphCreation()
          {
            auto world = net::MpiCommunicator::World();
            if(world.Size() >= 4)
            {
              auto graph = world.Graph({{1}, {0, 2, 3}, {1, 3}, {1, 2}});
              auto const neighbors = graph.GetNeighbors();
              if(graph.Rank() == 0)
              {
                CPPUNIT_ASSERT_EQUAL(size_t(1), neighbors.size());
                CPPUNIT_ASSERT_EQUAL(int(1), neighbors.front());
              }
              else if(graph.Rank() == 1)
              {
                CPPUNIT_ASSERT_EQUAL(size_t(3), neighbors.size());
                CPPUNIT_ASSERT_EQUAL(int(0), neighbors[0]);
                CPPUNIT_ASSERT_EQUAL(int(2), neighbors[1]);
                CPPUNIT_ASSERT_EQUAL(int(3), neighbors[2]);
              }
              else if(graph.Rank() == 2)
              {
                CPPUNIT_ASSERT_EQUAL(size_t(2), neighbors.size());
                CPPUNIT_ASSERT_EQUAL(int(1), neighbors[0]);
                CPPUNIT_ASSERT_EQUAL(int(3), neighbors[1]);
              }
              else if(graph.Rank() == 3)
              {
                CPPUNIT_ASSERT_EQUAL(size_t(2), neighbors.size());
                CPPUNIT_ASSERT_EQUAL(int(1), neighbors[0]);
                CPPUNIT_ASSERT_EQUAL(int(2), neighbors[1]);
              }
              else
              {
                CPPUNIT_ASSERT_EQUAL(size_t(0), neighbors.size());
              }
            }
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (MPINodeParallelizationTests);
    }
  }
}

#endif  // ONCE
