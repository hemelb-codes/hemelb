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
      using namespace hemelb::redblood;
      class BufferTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (BufferTests);
            CPPUNIT_TEST (testDrop);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            cylinder.normal = LatticePosition(1, 0, 0);
            cylinder.origin = LatticePosition(2, 2, 2);
            cylinder.radius = 1.5;
            buffer = std::make_shared<buffer::Buffer>(cylinder);

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
            buffer->SetOffset(LatticePosition(-10, 0, 10));

            // justDropped should be invalid since nothing was dropped
            CPPUNIT_ASSERT(not buffer->GetJustDropped());

            // Now make sure drop statements yield expected result
            for(auto i = 0; i < cells.size(); ++i)
            {
              auto cell = buffer->drop();
              auto expected_pos = positions[i] + buffer->GetOffset() - cylinder.origin;
              CPPUNIT_ASSERT(cell == cells[i]);
              CPPUNIT_ASSERT(is_zero(cell->GetBarycenter() - expected_pos));
              CPPUNIT_ASSERT(buffer->GetJustDropped() == cell);
            }

            // drop should throw if there are no cells.
            CPPUNIT_ASSERT_THROW(buffer->drop(), Exception);
          }

        private:
          std::shared_ptr<buffer::Buffer> buffer;
          std::vector<std::shared_ptr<Cell> > cells;
          Cylinder cylinder;
      };


      CPPUNIT_TEST_SUITE_REGISTRATION(BufferTests);
    }
  }
}

#endif  // ONCE
