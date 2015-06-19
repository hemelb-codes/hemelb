//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_VERTEXBAG_H
#define HEMELB_UNITTESTS_REDBLOOD_VERTEXBAG_H

#include <cppunit/TestFixture.h>
#include "redblood/VertexBag.h"

#include "redblood/VertexBag.cc"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class VertexBagTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (VertexBagTests);
          CPPUNIT_TEST (testConstruction);
          CPPUNIT_TEST (testProcsAffectedBySingleProc);
          CPPUNIT_TEST (testProcsAffectedByMultiProcs);
          CPPUNIT_TEST (testSplittingCellSingleProc);
          CPPUNIT_TEST (testSplittingCellMultiProcs);CPPUNIT_TEST_SUITE_END();

        public:
          void testConstruction()
          {
            auto const cell = std::make_shared<Cell>(icoSphere());
            VertexBag bag(cell);
            CPPUNIT_ASSERT(bag.GetTag() == cell->GetTag());
            CPPUNIT_ASSERT_EQUAL(size_t(0), bag.GetVertices().size());
          }

          void testProcsAffectedBySingleProc()
          {
            LatticePosition const X(10, 10, 10);
            auto const high_quadrant = procsAffectedByPosition(VertexBagTests::proc_at_pos, X);
            CPPUNIT_ASSERT_EQUAL(size_t(1), high_quadrant.size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), high_quadrant.count(0));
            auto const low_quadrant = procsAffectedByPosition(VertexBagTests::proc_at_pos, -X);
            CPPUNIT_ASSERT_EQUAL(size_t(1), low_quadrant.size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), low_quadrant.count(7));

            auto const avoid = procsAffectedByPosition(VertexBagTests::proc_at_pos, -X, 7);
            CPPUNIT_ASSERT_EQUAL(size_t(0), avoid.size());
          }

          void testProcsAffectedByMultiProcs()
          {
            LatticePosition const X(10, 0, 10);
            auto const two_quadrants = procsAffectedByPosition(VertexBagTests::proc_at_pos, X);
            CPPUNIT_ASSERT_EQUAL(size_t(2), two_quadrants.size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), two_quadrants.count(0));
            CPPUNIT_ASSERT_EQUAL(size_t(1), two_quadrants.count(2));

            LatticePosition const Y(0, 0, 10);
            auto const four_quadrants = procsAffectedByPosition(VertexBagTests::proc_at_pos, Y);
            CPPUNIT_ASSERT_EQUAL(size_t(4), four_quadrants.size());
            for (proc_t i(0); i < proc_t(4); ++i)
            {
              CPPUNIT_ASSERT_EQUAL(size_t(1), four_quadrants.count(i));
            }

            LatticePosition const Z(0, 0, 0);
            auto const all_quadrants = procsAffectedByPosition(VertexBagTests::proc_at_pos, Z);
            CPPUNIT_ASSERT_EQUAL(size_t(8), all_quadrants.size());
            for (proc_t i(0); i < proc_t(8); ++i)
            {
              CPPUNIT_ASSERT_EQUAL(size_t(1), all_quadrants.count(i));
            }
          }

          void testSplittingCellSingleProc()
          {
            auto const cell = std::make_shared<Cell>(icoSphere());
            *cell *= 5e0;
            *cell += LatticePosition(50, 50, 50) - cell->GetBarycenter();
            auto splits = splitVertices(VertexBagTests::proc_at_pos, cell);
            CPPUNIT_ASSERT_EQUAL(size_t(1), splits.size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), splits.count(0));
            auto const &expected = cell->GetVertices();
            auto const &actual = splits[0]->GetVertices();
            CPPUNIT_ASSERT_EQUAL(expected.size(), actual.size());
            for (size_t i(0); i < expected.size(); ++i)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i].x, actual[i].x, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i].y, actual[i].y, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[i].z, actual[i].z, 1e-8);
            }
          }

          void testSplittingCellMultiProcs()
          {
            auto const cell = std::make_shared<Cell>(icoSphere());
            *cell *= 5e0;
            *cell += LatticePosition(1, 0, 1) * 50e0 - cell->GetBarycenter();
            auto splits = splitVertices(VertexBagTests::proc_at_pos, cell);
            CPPUNIT_ASSERT_EQUAL(size_t(2), splits.size());
            CPPUNIT_ASSERT_EQUAL(size_t(1), splits.count(0));
            CPPUNIT_ASSERT_EQUAL(size_t(1), splits.count(2));
            // Four points are in common since their y coordinates is zero
            CPPUNIT_ASSERT_EQUAL(cell->GetVertices().size() + 4,
                                 splits[0]->GetVertices().size() + splits[2]->GetVertices().size());
          }

          static proc_t proc_at_pos(LatticePosition const &position)
          {
            return (position.x > 0 ?
              0 :
              1) + (position.y > 0 ?
              0 :
              1) * 2 + (position.z > 0 ?
              0 :
              1) * 4;
          }
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (VertexBagTests);
    }
  }
}

#endif  // ONCE
