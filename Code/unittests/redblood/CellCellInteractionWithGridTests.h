//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLCELL_INTERACTION_WITH_GRID_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLCELL_INTERACTION_WITH_GRID_TESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellCellInteractionWithGridTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE (CellCellInteractionWithGridTests);
          CPPUNIT_TEST (testInteraction);CPPUNIT_TEST_SUITE_END();

          PhysicalDistance const cutoff = 5.0;
          PhysicalDistance const halo = 2.0;

        public:
          void testInteraction();

        private:
          virtual size_t CubeSize() const
          {
            return 32 + 2;
          }
      };

      void CellCellInteractionWithGridTests::testInteraction()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Place two nodes close enough for interactions
        LatticePosition const n0(15 - 0.1, 15.5, 15.5);
        LatticePosition const n1(15 + 0.1, 15.5, 15.5);
        (*cells.begin())->GetVertices().front() = n0;
        (*std::next(cells.begin()))->GetVertices().front() = n1;

        // Set forces to zero
        helpers::ZeroOutFOld(latDat);

        // Finds pairs, computes interaction, spread forces to lattice
        addCell2CellInteractions(DivideConquerCells(cells, cutoff, halo),
                                 Node2NodeForce(1.0, halo),
                                 stencil::types::FOUR_POINT,
                                 *latDat);

        // By symmetry, there are no forces on the lattice points equidistant from
        // the nodes
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 14, 14).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 16, 16).GetForce()));

        // There are non-zero opposite forces on the following nodes
        CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(14, 15, 15).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(16, 15, 15).GetForce()
            + latDat->GetSite(14, 15, 15).GetForce()));
        // The forces at (14, 15, 15) should be  in direction (-1, 0, 0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(latDat->GetSite(14, 15, 15).GetForce().Dot(LatticePosition(-1,
                                                                                                0,
                                                                                                0)),
                                     std::abs(latDat->GetSite(14, 15, 15).GetForce().x),
                                     1e-8);

        // There are non-zero opposite forces on the following nodes
        CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(13, 14, 14).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(17, 14, 14).GetForce()
            + latDat->GetSite(13, 14, 14).GetForce()));
        CPPUNIT_ASSERT_DOUBLES_EQUAL(latDat->GetSite(13, 14, 14).GetForce().Dot(LatticePosition(-1,
                                                                                                0,
                                                                                                0)),
                                     std::abs(latDat->GetSite(13, 14, 14).GetForce().x),
                                     1e-8);

        // This node is too far away
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(12, 15, 15).GetForce()));
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellCellInteractionWithGridTests);
    }
  }
}

#endif  // ONCE
