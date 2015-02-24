//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_PARTICLESHUFFLER_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_PARTICLESHUFFLER_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/parallelization/ParticleShuffler.h"
// to unittest packing and unpacking in anonymous namespace
#include "redblood/parallelization/ParticleShuffler.cc"
#include <algorithm>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      //! Mock cell for ease of use
      class Shuffler : public ParticleShuffler
      {
        public:
          Shuffler(proc_t arank)
            : ParticleShuffler({arank == 0 ? 1: 0}), rank(arank)
          {
          }

          //! Figures out who owns an object
          std::tuple<bool, proc_t> GetCellOwner(CellContainer::const_reference cell) const override
          {
            auto const in0 = cell->GetBarycenter().GetMagnitude() < 5.0;
            return {not ((rank == 0) xor in0), in0 ? 0: 1};
          }
          //! Const reference to cells owned by this proc
          CellContainer const& GetOwnedCells() const override
          {
            return cells;
          }
          //! Returns an empty cell that will be filled by unpacking a message from another proc
          CellContainer::value_type GetEmptyCell() const override
          {
            CellContainer::value_type result(new Cell(*refcell));
            result->GetVertices().clear();
            result->SetScale(2e-12);
            return result;
          }
          //! Adds a cell as owned
          virtual void AddToOwnedCells(CellContainer::value_type&& cell) override
          {
            cells.emplace(std::move(cell));
          }
          //! Removes owned cell
          virtual void RemoveFromOwnedCells(CellContainer::const_reference cell) override
          {
            auto found = cells.find(cell);
            assert(found != cells.end());
            cells.erase(found);
          }

          proc_t rank = 0;
          CellContainer cells;
          CellContainer::value_type refcell;
      };

      class ParticleShufflerTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (ParticleShufflerTests);
          CPPUNIT_TEST (testPackUnpackCell);
          CPPUNIT_TEST (testNoMovement);
          CPPUNIT_TEST (testShuffleOneCell);
          CPPUNIT_TEST (testShuffleTwoCellsAndBack);
          CPPUNIT_TEST (testSloshing);
          CPPUNIT_TEST_SUITE_END();

        public:
          ParticleShufflerTests() : shuffler0(0), shuffler1(1)
          {
          }

          void setUp()
          {
            // adds two cells to shuffler0 only
            auto templateMesh = icoSphere();
            shuffler0.cells.emplace(new Cell(templateMesh.GetVertices(), templateMesh));
            templateMesh += 2.5;
            shuffler0.cells.emplace(new Cell(templateMesh.GetVertices(), templateMesh));

            shuffler0.refcell = *shuffler0.cells.begin();
            shuffler1.refcell = *shuffler0.cells.begin();

            messageSizeFor0 = 0;
            messageSizeFor1 = 0;
            nextMessageSizeFor0 = 0;
            nextMessageSizeFor1 = 0;
            messageFor0.resize(1);
            messageFor1.resize(1);
          }


          void testPackUnpackCell();
          void testNoMovement();
          void testShuffleOneCell();
          void testShuffleTwoCellsAndBack();
          void testSloshing();

        private:
          enum class Stepper : unsigned
          {
            IDENTIFY,
            NEXTCELLMESSAGE,
            SETCELLMESSAGE,
            PACK, // should be send in one phase
            UNPACK, // should be receive in same phase
            NEXT,  // moves to next phase/LB step
            END = NEXT,
            NONE,
            START = IDENTIFY
          };
          // Runs through the cycle and check no movement is made
          void checkNoMovement();
          void stepper(Stepper const &start, Stepper const & last = Stepper::NONE);

          Shuffler shuffler0;
          Shuffler shuffler1;
          size_t messageSizeFor0 = 0;
          size_t messageSizeFor1 = 0;
          size_t nextMessageSizeFor0 = 0;
          size_t nextMessageSizeFor1 = 0;
          std::vector<int8_t> messageFor0;
          std::vector<int8_t> messageFor1;
      };

      void ParticleShufflerTests :: testPackUnpackCell()
      {
        // Cell to pack and later unpack
        auto const cell = *shuffler0.cells.begin();
        // First figure out size
        auto const bufferSize = cellPackSize(cell);
        CPPUNIT_ASSERT(bufferSize > 0);
        // Create a larger buffer, just in case
        std::vector<int8_t> buffer(bufferSize * 2, 0);
        // pack the cell
        int8_t* endBuffer = (buffer.data() << cell);
        decltype(bufferSize) const diff(endBuffer - buffer.data());
        CPPUNIT_ASSERT_EQUAL(diff, bufferSize);

        CellContainer::value_type newCell(new Cell(*cell));
        newCell->GetVertices().clear();
        newCell->SetScale(2e12);
        CPPUNIT_ASSERT(newCell->GetNumberOfNodes() != cell->GetNumberOfNodes());
        CPPUNIT_ASSERT(std::abs(newCell->GetScale() - cell->GetScale()) > 0.1);
        endBuffer = buffer.data();
        endBuffer = (buffer.data() >> newCell);
        CPPUNIT_ASSERT_EQUAL(diff, bufferSize);
        CPPUNIT_ASSERT_EQUAL(newCell->GetNumberOfNodes(), cell->GetNumberOfNodes());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(newCell->GetScale(), cell->GetScale(), 1e-8);
        for(size_t i(0); i < newCell->GetNumberOfNodes(); ++i) {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(newCell->GetVertices()[i].x, cell->GetVertices()[i].x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(newCell->GetVertices()[i].y, cell->GetVertices()[i].y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(newCell->GetVertices()[i].z, cell->GetVertices()[i].z, 1e-8);
        }
      }

      void ParticleShufflerTests :: stepper(Stepper const & start, Stepper const & last)
      {
        auto doThisStep = [&start, &last](Stepper const &step)
        {
          return last == Stepper::NONE ? start == step: (start <= step and last >= step);
        };
        if(doThisStep(Stepper::IDENTIFY))
        {
          shuffler0.IdentifyOutOfBounds();
          shuffler1.IdentifyOutOfBounds();
        }
        if(doThisStep(Stepper::NEXTCELLMESSAGE))
        {
          nextMessageSizeFor1 = shuffler0.GetNextCellMessageSize(1);
          nextMessageSizeFor0 = shuffler1.GetNextCellMessageSize(0);
        }
        if(doThisStep(Stepper::SETCELLMESSAGE))
        {
          shuffler0.SetThisCellMessageSize(1, messageSizeFor0);
          shuffler1.SetThisCellMessageSize(0, messageSizeFor1);
        }
        if(doThisStep(Stepper::PACK))
        {
          messageFor0.resize(std::max(1ul, messageSizeFor0), 0);
          messageFor1.resize(std::max(1ul, messageSizeFor1), 0);
          int8_t* const advanced1 = shuffler0.Pack(1, messageFor1.data());
          int8_t* const advanced0 = shuffler1.Pack(0, messageFor0.data());
          CPPUNIT_ASSERT_EQUAL(int(messageSizeFor0), int(advanced0 - messageFor0.data()));
          CPPUNIT_ASSERT_EQUAL(int(messageSizeFor1), int(advanced1 - messageFor1.data()));
        }
        if(doThisStep(Stepper::UNPACK))
        {
          CPPUNIT_ASSERT_EQUAL(messageFor0.size(), std::max(1ul, messageSizeFor0));
          CPPUNIT_ASSERT_EQUAL(messageFor1.size(), std::max(1ul, messageSizeFor1));
          int8_t* const advanced0 = shuffler0.Unpack(1, messageFor0.data());
          int8_t* const advanced1 = shuffler1.Unpack(0, messageFor1.data());
          CPPUNIT_ASSERT_EQUAL(int(advanced0 - messageFor0.data()), int(messageSizeFor0));
          CPPUNIT_ASSERT_EQUAL(int(advanced1 - messageFor1.data()), int(messageSizeFor1));
        }
        if(doThisStep(Stepper::NEXT))
        {
          shuffler0.Next();
          shuffler1.Next();
          messageFor0.clear();
          messageFor1.clear();
          std::swap(nextMessageSizeFor0, messageSizeFor0);
          std::swap(nextMessageSizeFor1, messageSizeFor1);
        }
      }

      void ParticleShufflerTests :: checkNoMovement()
      {
        auto const cellSize0 = shuffler0.cells.size();
        auto const cellSize1 = shuffler1.cells.size();

        stepper(Stepper::IDENTIFY);
        stepper(Stepper::NEXTCELLMESSAGE);
        CPPUNIT_ASSERT_EQUAL(nextMessageSizeFor0, size_t(0));
        CPPUNIT_ASSERT_EQUAL(nextMessageSizeFor1, size_t(0));

        stepper(Stepper::SETCELLMESSAGE, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);

        // Ok now on second step
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);
      }

      void ParticleShufflerTests :: testNoMovement()
      {
        // Runs through the cycle twice, ensuring that the process is repeatable
        checkNoMovement();
        checkNoMovement();
      }

      void ParticleShufflerTests:: testShuffleOneCell()
      {
        auto const cellSize0 = shuffler0.cells.size();
        auto const cellSize1 = shuffler1.cells.size();

        checkNoMovement();

        // Moves one particle
        **shuffler0.cells.begin() += LatticePosition(50, 50, 50);
        size_t const messageSize = cellPackSize(*shuffler0.cells.begin());
        stepper(Stepper::IDENTIFY);
        stepper(Stepper::NEXTCELLMESSAGE);
        CPPUNIT_ASSERT(nextMessageSizeFor1 != 0);
        // One cell will be moved in *next* iteration ...
        CPPUNIT_ASSERT_EQUAL(nextMessageSizeFor1, messageSize);
        CPPUNIT_ASSERT_EQUAL(nextMessageSizeFor0, size_t(0));

        // ... where *next* iteration means not this one
        stepper(Stepper::SETCELLMESSAGE, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);

        // Allright, this is the iteration where stuff happens!
        CPPUNIT_ASSERT_EQUAL(messageSizeFor1, messageSize);
        CPPUNIT_ASSERT_EQUAL(messageSizeFor0, size_t(0));
        stepper(Stepper::IDENTIFY, Stepper::PACK);
        CPPUNIT_ASSERT_EQUAL(messageFor1.size(), messageSize);
        CPPUNIT_ASSERT_EQUAL(messageFor0.size(), size_t(1));
        // unpack should copy cell to one (but cell still in zero)
        stepper(Stepper::UNPACK);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);
        // next should remove cell from zero
        stepper(Stepper::NEXT);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 1);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);

        // Checking nothing further happens
        checkNoMovement();
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 1);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);
      }

      void ParticleShufflerTests:: testShuffleTwoCellsAndBack()
      {
        auto const cellSize0 = shuffler0.cells.size();
        auto const cellSize1 = shuffler1.cells.size();

        checkNoMovement();

        // Move two particles at a time
        **shuffler0.cells.begin() += LatticePosition(50, 50, 50);
        **std::next(shuffler0.cells.begin()) += LatticePosition(50, 50, 50);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 2);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 2);

        checkNoMovement();

        // Move one particles at a time
        **shuffler1.cells.begin() -= LatticePosition(50, 50, 50);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 2);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 2);
        **std::next(shuffler1.cells.begin()) -= LatticePosition(50, 50, 50);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 1);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);

        checkNoMovement();
      }

      void ParticleShufflerTests:: testSloshing()
      {
        auto const cellSize0 = shuffler0.cells.size();
        auto const cellSize1 = shuffler1.cells.size();

        checkNoMovement();

        // Particle moves back to shuffler0 before transit to shuffler1 is complete
        // However, the transit must still take place since communication channels are prepped for
        // it.
        **shuffler0.cells.begin() += LatticePosition(50, 50, 50);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);
        **shuffler0.cells.begin() -= LatticePosition(50, 50, 50);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 1);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);

        // Now the particle will actually move back
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0 - 1);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1 + 1);
        stepper(Stepper::START, Stepper::END);
        CPPUNIT_ASSERT_EQUAL(shuffler0.cells.size(), cellSize0);
        CPPUNIT_ASSERT_EQUAL(shuffler1.cells.size(), cellSize1);

        checkNoMovement();
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (ParticleShufflerTests);
    }
  }
}

#endif  // ONCE
