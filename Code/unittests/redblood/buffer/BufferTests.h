//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_BUFFER_BUFFER_H
#define HEMELB_UNITTESTS_REDBLOOD_BUFFER_BUFFER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/buffer/Buffer.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class OpenBuffer : public buffer::Buffer
      {
        public:
          OpenBuffer(std::shared_ptr<Cylinder> cyl, CellContainer const& cells = CellContainer())
            : buffer::Buffer(cyl, cells)
          {
          }
          OpenBuffer(Cylinder const & cyl, CellContainer const& cells = CellContainer())
            : buffer::Buffer(cyl, cells)
          {
          }

          CellContainer::value_type GetJustDropped() const
          {
            return justDropped;
          }
          LatticeDistance GetOffset() const
          {
            return offset;
          }
          LatticeDistance GetLastZ() const
          {
            return lastZ;
          }
          void SetOffset(LatticeDistance const &off)
          {
            offset = off;
          }
          void SetInteraction(LatticeDistance const &off)
          {
            interactionRadius = off;
          }
      };

      using namespace hemelb::redblood;
      class BufferTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (BufferTests);
            CPPUNIT_TEST (testDrop);
            CPPUNIT_TEST (testOffsetUpdateNoCell);
            CPPUNIT_TEST (testOffsetUpdateNoDropped);
            CPPUNIT_TEST (testOffsetUpdateDroppedMovedPerpendicular);
            CPPUNIT_TEST (testOffsetUpdateDroppedMovedBackward);
            CPPUNIT_TEST (testOffsetUpdateDroppedMovedForward);
            CPPUNIT_TEST (testOutsideInteractionRange);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            cylinder.normal = LatticePosition(1, 0, 0);
            cylinder.origin = LatticePosition(2, 2, 2);
            cylinder.radius = 1.5;
            buffer = std::make_shared<OpenBuffer>(cylinder);

            auto cell = std::make_shared<Cell>(icoSphere());
            *cell -= cell->GetBarycenter();
            cells.push_back(std::make_shared<Cell>(*cell));
            buffer->insert(cells.back());
            *cell += cylinder.normal * 2.0;
            cells.push_back(std::make_shared<Cell>(*cell));
            buffer->insert(cells.back());
            *cell += cylinder.normal * 4.0;
            cells.push_back(std::make_shared<Cell>(*cell));
            buffer->insert(cells.back());
          }

          void tearDown()
          {
          }

          // Drops cell closest to drop point and translates it to LB coords
          void testDrop()
          {
            // Gets original positions
            // We know the cells are in cells are sorted
            auto get_pos = [](CellContainer::value_type a)
            {
              return a->GetBarycenter();
            };
            std::vector<LatticePosition> positions(cells.size());
            std::transform(cells.begin(), cells.end(), positions.begin(), get_pos);

            // Change offset to make it interesting
            buffer->SetOffset(10);

            // justDropped should be invalid since nothing was dropped
            CPPUNIT_ASSERT(not buffer->GetJustDropped());

            // Now make sure drop statements yield expected result
            for(auto i = 0; i < cells.size(); ++i)
            {
              auto cell = buffer->drop();
              auto expected_pos
                = positions[i] + cylinder.normal * buffer->GetOffset() + cylinder.origin;
              CPPUNIT_ASSERT(cell == cells[i]);
              CPPUNIT_ASSERT(is_zero(cell->GetBarycenter() - expected_pos));
              CPPUNIT_ASSERT(buffer->GetJustDropped() == cell);
            }

            // drop should throw if there are no cells.
            CPPUNIT_ASSERT_THROW(buffer->drop(), Exception);
          }

          void testOffsetUpdateNoCell()
          {
            // Remove all cells
            for(size_t i=0; i < cells.size(); ++i)
            {
              buffer->drop();
            }
            // Offset should be reset to zero
            buffer->SetOffset(10);
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, buffer->GetOffset(), 1e-8);
          }

          void testOffsetUpdateNoDropped()
          {
            // set first cell to something not zero, to make the test interesting
            *cells.front() += cylinder.normal * 1.0;
            buffer->updateOffset();
            // updating should set offet such that dropping the cell will place it at the origin
            auto const cell = buffer->drop();
            auto const barycenter = cell->GetBarycenter();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cylinder.origin.x, barycenter.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cylinder.origin.y, barycenter.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cylinder.origin.z, barycenter.z, 1e-8);
          }

          // Moving perpendicular to the normal should have no effect
          void testOffsetUpdateDroppedMovedPerpendicular()
          {
            auto cell = buffer->drop();
            buffer->updateOffset();
            auto const startOffset = buffer->GetOffset();

            *cell += LatticePosition(0, 1, 0) * 1.50;
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset, buffer->GetOffset(), 1e-8);

            *cell += LatticePosition(0, 0, 1) * 15.0;
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset, buffer->GetOffset(), 1e-8);
          }

          // Moving backwards when cell are interacting should keep offset fixed (buffer is
          // immoveable)
          void testOffsetUpdateDroppedMovedBackward()
          {
            // Special case only if dropped and next cell interact
            buffer->SetInteraction(15e0);
            buffer->updateOffset();
            auto cell = buffer->drop();
            auto const startOffset = buffer->GetOffset();

            *cell += cylinder.normal * 15.0;
            buffer->updateOffset();
            auto const expected = buffer->GetOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset, expected, 1e-8);
          }

          // When interacting, cannot move by more than what the cells moves forward
          // And cannot move next cell beyond drop point
          void testOffsetUpdateDroppedMovedForward()
          {
            // Special case only if dropped and next cell interact
            buffer->SetInteraction(15e0);
            buffer->updateOffset();
            auto cell = buffer->drop();
            auto const startOffset = buffer->GetOffset();

            auto const delta = 0.1;
            *cell -= cylinder.normal * delta;
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset-delta, buffer->GetOffset(), 1e-8);

            *cell -= cylinder.normal * 2.0;
            buffer->updateOffset();
            auto const nextDrop = buffer->drop();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(cylinder.origin.x, nextDrop->GetBarycenter().x, 1e-8);
          }

          void testOutsideInteractionRange()
          {
            // Special case only if dropped and next cell interact
            buffer->SetInteraction(1e0);

            // Arrange cell positions first
            *(cells[0]) += cylinder.normal *  0.0 - cells[0]->GetBarycenter();
            *(cells[1]) += cylinder.normal *  5.0 - cells[1]->GetBarycenter();
            *(cells[2]) += cylinder.normal * 10.0 - cells[2]->GetBarycenter();

            buffer->updateOffset();
            auto cell = buffer->drop();
            auto const startOffset = buffer->GetOffset();

            CPPUNIT_ASSERT(cell == cells[0]);

            // Move forward the first cell, but not so little that they would interact
            *cell -= cylinder.normal * 10.0;
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset - 5e0, buffer->GetOffset(), 1e-8);

            // Now drop another cell. Offset is limited by interactions.
            cell = buffer->drop();
            CPPUNIT_ASSERT(cell == cells[1]);
            buffer->updateOffset();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset - 5e0 - 4e0, buffer->GetOffset(), 1e-8);
          }

        private:
          std::shared_ptr<OpenBuffer> buffer;
          std::vector<std::shared_ptr<Cell> > cells;
          Cylinder cylinder;
      };



      CPPUNIT_TEST_SUITE_REGISTRATION(BufferTests);
    }
  }
}

#endif  // ONCE
