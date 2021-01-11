// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPINODEPARALLELIZATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_MPINODEPARALLELIZATIONTESTS_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "net/MpiCommunicator.h"
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
          CPPUNIT_TEST (testNeighborhoodGraphCreation);CPPUNIT_TEST_SUITE_END();

        public:
          //! Test creation of mpi graph topology
          //! \details The test topology includes only the first four nodes:
          //!
          //!            2
          //!           / \
          //!          /   \
          //!  0 ---- 1 ---- 3  4  5  6 ...
          //!
          //! This topology checks we have nodes with different numbers of neighbors, cycles, and
          //! empty neighborhoods.
          void testNeighborhoodGraphCreation();
      };

      void MPINodeParallelizationTests::testNeighborhoodGraphCreation()
      {
        auto world = net::MpiCommunicator::World();
        if (world.Size() >= 4)
        {
          std::vector<std::vector<int>> vertices { { 1 }, { 0, 2, 3 }, { 1, 3 }, { 1, 2 } };
          for (int i(4); i < world.Size(); ++i)
          {
            vertices.push_back(std::vector<int> { });
          }
          auto graph = world.Graph(vertices);
          auto const neighbors = graph.GetNeighbors();
          CPPUNIT_ASSERT_EQUAL(vertices[graph.Rank()].size(), neighbors.size());
          for (auto const item : util::zip(vertices[graph.Rank()], neighbors))
          {
            CPPUNIT_ASSERT_EQUAL(std::get<0>(item), std::get<1>(item) + 1);
          }
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (MPINodeParallelizationTests);
    }
  }
}

#endif  // ONCE
