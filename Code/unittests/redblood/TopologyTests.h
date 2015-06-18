//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_TOPOLOGY_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_TOPOLOGY_TESTS_H

#include <sstream>
#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "resources/Resource.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class TopologyTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (TopologyTests);
          CPPUNIT_TEST (testNodeToVertex);
          CPPUNIT_TEST (testFacetNeighbors);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            std::string filename = resources::Resource("red_blood_cube.txt").Path();
            mesh = readMesh(filename);
            // Checks the mesh input makes sense
            CPPUNIT_ASSERT(mesh->vertices.size() == 8);
            CPPUNIT_ASSERT(mesh->facets.size() == 12);

            // Creates topology
            topo = std::shared_ptr<MeshTopology>(new MeshTopology(*mesh));
          }

          void tearDown()
          {
          }

          void testNodeToVertex()
          {
            CPPUNIT_ASSERT(topo->vertexToFacets.size() == 8);

            // expected[vertex] = {nfacets, facet indices}
            unsigned int expected[8][6] = { { 4, 0, 1, 6, 9 }, { 5, 0, 2, 3, 8, 9 }, { 4,
                                                                                       0,
                                                                                       1,
                                                                                       3,
                                                                                       10 },
                                            { 5, 1, 6, 7, 10, 11 }, { 5, 4, 6, 7, 8, 9 }, { 4,
                                                                                            2,
                                                                                            4,
                                                                                            5,
                                                                                            8 },
                                            { 5, 2, 3, 5, 10, 11 }, { 4, 4, 5, 7, 11 }, };

            for (unsigned vertex(0); vertex < 8; ++vertex)
            {
              std::set<size_t> facets = topo->vertexToFacets[vertex];
              CPPUNIT_ASSERT(facets.size() == expected[vertex][0]);

              for (size_t facet(1); facet <= expected[vertex][0]; ++facet)
              {
                CPPUNIT_ASSERT(facets.count(expected[vertex][facet]) == 1);
              }
            }
          }

          void testFacetNeighbors()
          {
            CPPUNIT_ASSERT(topo->facetNeighbors.size() == 12);

            // expected[facet] = {neighbor indices}
            size_t expected[12][3] = { { 1, 3, 9 }, { 0, 6, 10 }, { 3, 5, 8 }, { 0, 2, 10 }, { 5,
                                                                                               7,
                                                                                               8 },
                                       { 4, 2, 11 }, { 1, 7, 9 }, { 6, 4, 11 }, { 9, 4, 2 }, { 0,
                                                                                               6,
                                                                                               8 },
                                       { 1, 3, 11 }, { 10, 5, 7 }, };

            for (unsigned facet(0); facet < 12; ++facet)
            {
              std::array<size_t, 3> const &neighs = topo->facetNeighbors[facet];
              CPPUNIT_ASSERT(neighs.size() == 3);

              for (size_t neigh(0); neigh < 3; ++neigh)
              {
                CPPUNIT_ASSERT(contains(neighs, expected[facet][neigh]));
              }
            }
          }

          template<class T, unsigned long N>
          static bool contains(std::array<T, N> const &facets, T value)
          {
            for (unsigned i(0); i < N; ++i)
              if (facets[i] == value)
              {
                return true;
              }

            return false;
          }

        private:
          std::shared_ptr<MeshData> mesh;
          std::shared_ptr<MeshTopology> topo;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (TopologyTests);
    }
  }
}

#endif  // ONCE
