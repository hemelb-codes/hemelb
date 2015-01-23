//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_FACET_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_FACET_TESTS_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/Cell.impl.cc"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      // Tests functionality that is *not* part of the HemeLB API
      // Checks that we know how to compute geometric properties between facets
      // However, HemeLB only cares about energy and forces
      class FacetTests : public BasisFixture
      {
          CPPUNIT_TEST_SUITE (FacetTests);
          CPPUNIT_TEST (testNormal);
          CPPUNIT_TEST (testUnitNormal);
          CPPUNIT_TEST (testAngle);
          CPPUNIT_TEST (testCommonNodes);
          CPPUNIT_TEST (testOrientedAngle);
          CPPUNIT_TEST (testSingleNodes);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            BasisFixture::setUp();
            main.reset(new Facet(mesh, 0));
            neighbor.reset(new Facet(mesh, 3));
          }

          void tearDown()
          {
          }

          void testNormal()
          {
            {
              LatticePosition const actual = main->normal();
              LatticePosition const expected(0, 0, -1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
            {
              LatticePosition const actual = neighbor->normal();
              LatticePosition const expected(1, 1, 1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
          }

          void testUnitNormal()
          {
            {
              LatticePosition const actual = main->unitNormal();
              LatticePosition const expected(0, 0, -1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected));
            }
            {
              LatticePosition const actual = neighbor->unitNormal();
              LatticePosition const expected(1, 1, 1);
              CPPUNIT_ASSERT(helpers::is_zero(actual - expected.GetNormalised()));
            }
          }

          void testAngle()
          {
            Angle const actual0 = angle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual0 - std::acos(-1.0 / std::sqrt(3.))));

            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            Angle const actual1 = angle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 - 3.0 * PI / 4.0));
            mesh.vertices.back()[1] = 1e0;
          }

          void testCommonNodes()
          {
            IndexPair const nodes = commonNodes(*main, *neighbor);
            CPPUNIT_ASSERT(nodes.first == 1 and nodes.second == 2);
            LatticePosition const edge = commonEdge(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(edge - (mesh.vertices[1] - mesh.vertices[2])));
          }
          void testSingleNodes()
          {
            IndexPair const nodes = singleNodes(*main, *neighbor);
            // indices into the list of vertices of each facet,
            // not into list of all vertices.
            CPPUNIT_ASSERT(nodes.first == 0 and nodes.second == 1);
          }

          void testOrientedAngle()
          {
            Angle const actual0 = orientedAngle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual0 + std::acos(-1.0 / std::sqrt(3.))));

            // simpler angle
            mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
            Angle const actual1 = orientedAngle(*main, *neighbor);
            CPPUNIT_ASSERT(helpers::is_zero(actual1 + 3.0 * PI / 4.0));

            // change orientation <==> negative angle
            mesh.facets.front()[1] = 2;
            mesh.facets.front()[2] = 1;
            Angle const actual2 = orientedAngle(Facet(mesh, 0), Facet(mesh, 3));
            CPPUNIT_ASSERT(helpers::is_zero(actual2 + PI / 4.0));
            mesh.facets.front()[1] = 1;
            mesh.facets.front()[2] = 2;

            mesh.vertices.back()[1] = 1e0;
          }

        protected:
          LatticePosition nodes[4];
          std::set<size_t> main_indices;
          std::set<size_t> neighbor_indices;
          std::shared_ptr<Facet> main, neighbor;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (FacetTests);
    }
  }
}

#endif  // ONCE
