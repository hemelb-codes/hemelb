//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_BUFFER_COLUMNSTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_BUFFER_COLUMNSTESTS_H

#include "unittests/redblood/Fixtures.h"
#include <cppunit/TestFixture.h>

#include "redblood/buffer/Columns.h"
// includes functions from anonymous namespace
#include "redblood/buffer/Columns.cc"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class ColumnsTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (ColumnsTests);
          CPPUNIT_TEST (testIdentityRotationMatrix);
          CPPUNIT_TEST (testRotationMatrix90Degrees);
          CPPUNIT_TEST (testRotationMatrix45Degrees);
          CPPUNIT_TEST (testMaxExtensions);
          CPPUNIT_TEST (testIterator);
          CPPUNIT_TEST (testCellDrop);CPPUNIT_TEST_SUITE_END();

        public:
          ColumnsTests() :
              CppUnit::TestFixture(), colAxis(1, 0, 0), cellAxis(1, 1, 1), sep(1)
          {
          }

          void setUp()
          {
            cylinder = std::make_shared<Cylinder>();
            cylinder->normal = LatticePosition(0, 0, 1);
            cylinder->origin = LatticePosition(5, 5, 0);
            cylinder->radius = 10;
          }

          void testIdentityRotationMatrix()
          {
            using namespace hemelb::redblood::buffer;
            auto const r = rotMat(LatticePosition(1, 0, 0), LatticePosition(1, 0, 0));
            for (size_t i(0); i < 3; ++i)
              for (size_t j(0); j < 3; ++j)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(r[i][j], i == j ?
                  1e0 :
                  0e0,
                                             1e-8);
          }
          void testRotationMatrix90Degrees()
          {
            using namespace hemelb::redblood::buffer;
            auto const r = rotMat(LatticePosition(1, 0, 0), LatticePosition(0, 1, 0));
            // a0 maps to b0
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(1, 0, 0)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(1, 0, 0)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(1, 0, 0)).z, 1e-8);
            // a0.Cross(b0) maps to itsel   f
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, 0, 1)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, 0, 1)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(0, 0, 1)).z, 1e-8);
            // a0.Cross(a0.Cross(b0)) maps to b0.Cross(a0.Cross(b0))
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(0, -1, 0)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, -1, 0)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, -1, 0)).z, 1e-8);
          }
          void testRotationMatrix45Degrees()
          {
            using namespace hemelb::redblood::buffer;
            auto const r = rotMat(LatticePosition(1, 0, 0), LatticePosition(1, 1, 0));
            PhysicalDistance const sqrt2 = std::sqrt(2e0);
            // a0 maps to b0
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(sqrt2, 0, 0)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(sqrt2, 0, 0)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(sqrt2, 0, 0)).z, 1e-8);
            // a0.Cross(b0) maps to itself
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, 0, 1)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(0, 0, 1)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, (r * LatticePosition(0, 0, 1)).z, 1e-8);
            // a0.Cross(a0.Cross(b0)) maps to b0.Cross(a0.Cross(b0))
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(1, 1, 0)).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt2, (r * LatticePosition(1, 1, 0)).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0, (r * LatticePosition(1, 1, 0)).z, 1e-8);
          }

          void testMaxExtensions()
          {
            using namespace hemelb::redblood::buffer;
            auto const verts = tetrahedron().GetVertices();
            PhysicalDistance const s2 = std::sqrt(2e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(s2, maxExtension(verts, LatticePosition(1, 1, 0)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(s2, maxExtension(verts, LatticePosition(-1, -1, 0)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, maxExtension(verts, LatticePosition(1, 0, 0)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, maxExtension(verts, LatticePosition(0, 0, 1)), 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1, maxExtension(verts, LatticePosition(0, 0, 0.5)), 1e-8);
          }

          void testIterator()
          {
            using namespace hemelb::redblood::buffer;
            // create deformed mesh
            MeshData::Vertices verts = tetrahedron().GetVertices();
            verts.front() += LatticePosition(0.25, 0.25, 0.25);
            ColumnPositionIterator iterator(cylinder, verts, cellAxis, colAxis, sep);

            // Find first vector: should be colinear with column
            std::vector<LatticePosition> positions;
            positions.push_back(*iterator);
            CPPUNIT_ASSERT(is_in_cylinder(positions[0], verts));
            ++iterator;
            positions.push_back(*iterator);
            CPPUNIT_ASSERT(is_in_cylinder(*iterator, verts));
            const auto a0 = positions[1] - positions[0];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, a0.Cross(colAxis).GetMagnitude(), 1e-8);
            // going two back should be outside cylinder
            // Goin one back might not since checking size of mesh is only approximate
            CPPUNIT_ASSERT(not is_in_cylinder(positions[0] - a0 * 2.0, verts));

            // Now go forward until we change column
            do
            {
              ++iterator;
              positions.push_back(*iterator);
              CPPUNIT_ASSERT(is_in_cylinder(positions.back(), verts));
            }
            while ( (positions.back() - positions[0]).Cross(a0).GetMagnitude() < 1e-8);

            // a1 goes from last position in one column to first position in other
            auto const a1 = positions.back() - positions[positions.size() - 2];
            CPPUNIT_ASSERT(a1.Cross(a0).GetMagnitude() > 1e-8);
            CPPUNIT_ASSERT(a1.Cross(a0).Cross(cylinder->normal).GetMagnitude() < 1e-8);

            // Now look for first item outside current plane
            do
            {
              ++iterator;
              positions.push_back(*iterator);
              CPPUNIT_ASSERT(is_in_cylinder(positions.back(), verts));
            }
            while ( (positions.back() - positions[0]).Dot(cylinder->normal) < 1e-8);

            // Check next positions are same as previous but translated along normal
            // positions.back() == positions[i] + something * cylinder->normal.
            const auto N = positions.size();
            CPPUNIT_ASSERT(N > 1);
            for (size_t i = 0; i < 2 * N; ++i)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(positions[i].x, positions.back().x, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(positions[i].y, positions.back().y, 1e-8);
              ++iterator;
              positions.push_back(*iterator);
              CPPUNIT_ASSERT(is_in_cylinder(positions.back(), verts));
            }
          }

          void testCellDrop()
          {
            using namespace hemelb::redblood::buffer;
            auto templateCell = std::make_shared<Cell>(tetrahedron());
            templateCell->moduli.bending = 5.0;
            templateCell->moduli.volume = 3.0;
            auto const rotation = rotMat(cellAxis, colAxis);

            auto check = [templateCell, this, &rotation](std::shared_ptr<CellBase const> acell)
            {
              auto const cell = std::static_pointer_cast<Cell const, CellBase const>(acell);
              auto const b0 = templateCell->GetBarycenter();
              auto const b1 = cell->GetBarycenter();
              auto const& vertices0 = templateCell->GetVertices();
              auto const& vertices1 = cell->GetVertices();
              CPPUNIT_ASSERT_EQUAL(templateCell->GetNumberOfNodes(), cell->GetNumberOfNodes());
              for(site_t i(0); i < cell->GetNumberOfNodes(); ++i)
              {
                // Cells are rotated with respect to original input by a given value
                // They are translated by a value that changes at each iteration
                // The translation is not checked here.
                auto const expected = rotation * (vertices0[i] - b0);
                auto const actual = vertices1[i] - b1;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.x, actual.x, 1e-8);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.y, actual.y, 1e-8);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected.z, actual.z, 1e-8);
                CPPUNIT_ASSERT(this->is_in_cylinder(b1, cell->GetVertices()));
              }
              CPPUNIT_ASSERT(templateCell->GetTemplateMesh().isSameData(cell->GetTemplateMesh()));
              CPPUNIT_ASSERT_DOUBLES_EQUAL(
                  templateCell->moduli.bending, cell->moduli.bending, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(
                  templateCell->moduli.volume, cell->moduli.volume, 1e-8);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(
                  templateCell->moduli.surface, cell->moduli.surface, 1e-8);
            };
            ColumnCellDrop dropCell(cylinder, templateCell, cellAxis, colAxis, sep);
            for (size_t i(0); i < 50; ++i)
            {
              check(dropCell());
            }
          }
        private:
          std::shared_ptr<Cylinder> cylinder;
          LatticePosition const colAxis;
          LatticePosition const cellAxis;
          LatticeDistance const sep;

          bool is_in_cylinder(LatticePosition const &a, MeshData::Vertices const &verts) const
          {
            LatticePosition const barycenter = hemelb::redblood::barycenter(verts);
            LatticePosition const n0 = cylinder->normal.GetNormalised();
            if (a.Dot(cylinder->normal) < -1e-8)
            {
              return false;
            }
            for (auto const &v : verts)
            {
              LatticePosition const x = a + v - barycenter - cylinder->origin;
              if (x.Cross(n0).GetMagnitude() >= cylinder->radius)
              {
                return false;
              }
            }
            return true;
          }
          ;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (ColumnsTests);
    }
  }
}

#endif  // ONCE
