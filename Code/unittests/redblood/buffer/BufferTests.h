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
#include "unittests/helpers/Comparisons.h"
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
          OpenBuffer(std::shared_ptr<Cylinder> cyl, CellContainer const& cells = CellContainer()) :
              buffer::Buffer(cyl, cells)
          {
          }
          OpenBuffer(Cylinder const & cyl, CellContainer const& cells = CellContainer()) :
              buffer::Buffer(cyl, cells)
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
          double GetInteraction() const
          {
            return interactionRadius;
          }
          auto insert(CellContainer::value_type cell) -> decltype(CellContainer().insert(cell))
          {
            return buffer::Buffer::insert(cell);
          }
          void updateOffset()
          {
            return buffer::Buffer::updateOffset();
          }
          CellContainer::value_type drop()
          {
            return buffer::Buffer::drop();
          }
          CellContainer::value_type nearestCell() const
          {
            return buffer::Buffer::nearestCell();
          }
          CellContainer::value_type furthestCell() const
          {
            return buffer::Buffer::furthestCell();
          }
          bool isDroppablePosition(CellContainer::value_type const &candidate) const
          {
            return buffer::Buffer::isDroppablePosition(candidate);
          }
          size_t size() const
          {
            return virtuals.size();
          }
          void clear()
          {
            justDropped.reset();
            virtuals.clear();
          }
          void fillBuffer(size_t n)
          {
            buffer::Buffer::fillBuffer(n);
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
          CPPUNIT_TEST (testNearAndFar);
          CPPUNIT_TEST (testIsDroppable);
          CPPUNIT_TEST (testFillBuffer);
          CPPUNIT_TEST (testFillBufferThrowIfNoFunc);
          CPPUNIT_TEST (testNumberOfRequests);
          CPPUNIT_TEST (testDropsCellsWhenPossible);
          CPPUNIT_TEST (testNoDropsCellsWhenNotPossible);CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            cylinder.normal = LatticePosition(1, 0, 0);
            cylinder.origin = LatticePosition(2, 2, 2);
            cylinder.radius = 1.5;
            cylinder.length = 2.0;
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

            increment = cylinder.normal * 0.5;
            cell = cells[0]->clone();
            auto filler = [cell, this]()
            {
              currentFillerPosition += this->increment;
              CellContainer::value_type result(cell->clone().release());
              *result += currentFillerPosition;
              this->fillers.push_back(result);
              return result;
            };
            buffer->SetNewCellFunction(filler);
            buffer->SetInteraction(0.6);
            buffer->SetOffset(0);
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
            for (auto i = 0; i < cells.size(); ++i)
            {
              auto cell = buffer->drop();
              auto expected_pos = positions[i] + cylinder.normal * buffer->GetOffset()
                  + cylinder.origin;
              CPPUNIT_ASSERT(cell == cells[i]);
              CPPUNIT_ASSERT(helpers::is_zero(cell->GetBarycenter() - expected_pos));
              CPPUNIT_ASSERT(buffer->GetJustDropped() == cell);
            }

            // drop should throw if there are no cells.
            CPPUNIT_ASSERT_THROW(buffer->drop(), Exception);
          }

          void testOffsetUpdateNoCell()
          {
            // Remove all cells
            for (size_t i = 0; i < cells.size(); ++i)
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
            CPPUNIT_ASSERT_DOUBLES_EQUAL(startOffset - delta, buffer->GetOffset(), 1e-8);

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
            * (cells[0]) += cylinder.normal * 0.0 - cells[0]->GetBarycenter();
            * (cells[1]) += cylinder.normal * 5.0 - cells[1]->GetBarycenter();
            * (cells[2]) += cylinder.normal * 10.0 - cells[2]->GetBarycenter();

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

          void testNearAndFar()
          {
            auto const near = buffer->nearestCell();
            CPPUNIT_ASSERT(near.get() == cells.front().get());
            auto const far = buffer->furthestCell();
            CPPUNIT_ASSERT(far.get() == cells.back().get());
          }

          void testIsDroppable()
          {
            auto cell = buffer->nearestCell();
            *cell -= cell->GetBarycenter();
            // Center of cylinder should pass
            *cell += cylinder.normal * cylinder.length * 0.5;
            buffer->SetOffset(0);
            CPPUNIT_ASSERT(buffer->isDroppablePosition(cell));
            // Offset moves cylinder towards vascular system
            // Now cell is outside cylinder and outside vascular system
            buffer->SetOffset(cylinder.length);
            CPPUNIT_ASSERT(not buffer->isDroppablePosition(cell));
            // Move cell back into cylinder
            *cell -= cylinder.normal * cylinder.length;
            CPPUNIT_ASSERT(buffer->isDroppablePosition(cell));
            // Move cell down into vascular system
            *cell -= cylinder.normal * cylinder.length;
            CPPUNIT_ASSERT(buffer->isDroppablePosition(cell));
          }

          void testFillBufferThrowIfNoFunc()
          {
            auto buf = std::make_shared<OpenBuffer>(cylinder);
            CPPUNIT_ASSERT_THROW(buf->fillBuffer(0), Exception);
          }

          void testFillBuffer()
          {
            currentFillerPosition = LatticePosition(0, 0, 0);
            buffer->clear();
            CPPUNIT_ASSERT_EQUAL(0, static_cast<int>(buffer->size()));
            buffer->fillBuffer(0);
            // We should have enough to fill cylinder.length + interaction of these critters
            int const ncells = static_cast<int>( (std::ceil(cylinder.length
                + buffer->GetInteraction()) / 0.5));
            CPPUNIT_ASSERT_EQUAL(ncells, static_cast<int>(buffer->size()));
            // We can still add more if wanted. However, once full, only adds as many as requested.
            buffer->fillBuffer(0);
            CPPUNIT_ASSERT_EQUAL(ncells, static_cast<int>(buffer->size()));
            buffer->fillBuffer(2);
            CPPUNIT_ASSERT_EQUAL(ncells + 2, static_cast<int>(buffer->size()));
          }

          void testNumberOfRequests()
          {
            CPPUNIT_ASSERT_EQUAL(site_t(0), buffer->NumberOfRequests());
            buffer->requestNewCells(1);
            CPPUNIT_ASSERT_EQUAL(site_t(1), buffer->NumberOfRequests());
            buffer->requestNewCells(10);
            CPPUNIT_ASSERT_EQUAL(site_t(11), buffer->NumberOfRequests());
          }
          void testDropsCellsWhenPossible()
          {
            currentFillerPosition = LatticePosition(0, 0, 0);
            std::vector<CellContainer::value_type> transfered;
            auto addcell = [&transfered](CellContainer::value_type a)
            {
              transfered.push_back(a);
            };
            int const ncells = static_cast<int>( (std::ceil(cylinder.length
                + buffer->GetInteraction()) / 0.5));

            // There should be no drop if there are no requests
            buffer->clear();
            CPPUNIT_ASSERT(static_cast<int>(buffer->size()) == 0);
            (*buffer)(addcell);
            CPPUNIT_ASSERT_EQUAL(size_t(0), transfered.size());

            // Now are going to drop n cells
            buffer->requestNewCells(ncells - 1);
            (*buffer)(addcell);
            // We expect n cells along given direction were transfered
            CPPUNIT_ASSERT_EQUAL(size_t(ncells - 1), transfered.size());
            for (size_t i(0); i < transfered.size(); ++i)
            {
              // those cells should be those created when filling the buffer
              CPPUNIT_ASSERT(transfered[i].get() == fillers[i].get());
            }
            // we expect buffer should still be full
            CPPUNIT_ASSERT(static_cast<int>(buffer->size()) >= 1);
          }

          void testNoDropsCellsWhenNotPossible()
          {
            currentFillerPosition = LatticePosition(0, 0, 0);
            std::vector<CellContainer::value_type> transfered;
            auto addcell = [&transfered](CellContainer::value_type a)
            {
              transfered.push_back(a);
            };
            int const ncells = static_cast<int>( (std::ceil(cylinder.length
                + buffer->GetInteraction()) / 0.5));

            buffer->clear();
            buffer->requestNewCells(ncells);
            (*buffer)(addcell);
            // only dropped n - 1 because outside buffer
            CPPUNIT_ASSERT_EQUAL(size_t(ncells - 1), transfered.size());
            CPPUNIT_ASSERT_EQUAL(site_t(1), buffer->NumberOfRequests());
            CPPUNIT_ASSERT(static_cast<int>(buffer->size()) >= 1);
            // now move last cell to make space for another, but not enough
            * (transfered.back()) -= cylinder.normal * 0.500 * buffer->GetInteraction();
            (*buffer)(addcell);
            CPPUNIT_ASSERT_EQUAL(size_t(ncells - 1), transfered.size());
            CPPUNIT_ASSERT_EQUAL(site_t(1), buffer->NumberOfRequests());
            // now move last cell to make space for another
            * (transfered.back()) -= cylinder.normal * 0.501 * buffer->GetInteraction();
            buffer->requestNewCells(1);
            (*buffer)(addcell);
            CPPUNIT_ASSERT_EQUAL(size_t(ncells), transfered.size());
            CPPUNIT_ASSERT_EQUAL(site_t(1), buffer->NumberOfRequests());
            CPPUNIT_ASSERT(static_cast<int>(buffer->size()) >= 1);
          }

        private:
          std::shared_ptr<OpenBuffer> buffer;
          std::vector<std::shared_ptr<Cell> > cells;
          std::vector<CellContainer::value_type> fillers;
          LatticePosition currentFillerPosition;
          LatticePosition increment;
          Cylinder cylinder;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (BufferTests);
    }
  }
}

#endif  // ONCE
