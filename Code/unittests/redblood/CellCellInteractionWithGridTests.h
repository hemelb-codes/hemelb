// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLCELLINTERACTIONWITHGRIDTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLCELLINTERACTIONWITHGRIDTESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellCellInteractionWithGridTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE (CellCellInteractionWithGridTests);
          CPPUNIT_TEST (testInteraction<stencil::FourPoint> );
          CPPUNIT_TEST (testInteraction<stencil::CosineApprox> );
          CPPUNIT_TEST (testInteraction<stencil::ThreePoint> );
          CPPUNIT_TEST (testInteraction<stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

          LatticeDistance const cutoff = 5.0;
          LatticeDistance const halo = 2.0;

        public:
          template<class STENCIL> void testInteraction();

        private:
          virtual size_t CubeSize() const
          {
            return 32 + 2;
          }
      };

      template<class STENCIL> void CellCellInteractionWithGridTests::testInteraction()
      {
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Place two nodes close enough for interactions
        LatticePosition const n0 { 15 - 0.1, 15.5, 15.5 };
        LatticePosition const n1 { 15 + 0.1, 15.5, 15.5 };
        (*cells.begin())->GetVertices().front() = n0;
        (*std::next(cells.begin()))->GetVertices().front() = n1;

        // Set forces to zero
        helpers::ZeroOutFOld(latDat);

        // Finds pairs, computes interaction, spread forces to lattice
        addCell2CellInteractions<STENCIL>(DivideConquerCells(cells, cutoff, halo),
                                          Node2NodeForce(1.0, halo),
                                          *latDat);

        // By symmetry, there are no forces on the lattice points equidistant from
        // the nodes
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 14, 14).GetForce()));
        CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(15, 16, 16).GetForce()));

        // There are non-zero opposite forces on the following nodes
        size_t delta(1);
        for (; 2 * delta <= STENCIL::GetRange(); ++delta)
        {
          CPPUNIT_ASSERT(latDat->GetSite(14, 15, 15).GetForce().GetMagnitudeSquared() > 1e-8);
          auto const left = latDat->GetSite(15 + delta, 15, 15).GetForce();
          auto const right = latDat->GetSite(15 - delta, 15, 15).GetForce();
          CPPUNIT_ASSERT_DOUBLES_EQUAL(left.x, -right.x, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(left.y, -right.y, 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(left.z, -right.z, 1e-8);
          // The forces at (14, 15, 15) should be  in direction (-1, 0, 0)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(right.Dot( util::Vector3D<double>{-1, 0, 0}), std::abs(right.x), 1e-8);
        }
        // This node is too far away
        CPPUNIT_ASSERT(latDat->GetSite(15 + delta, 15, 15).GetForce().GetMagnitudeSquared() < 1e-8);
        CPPUNIT_ASSERT(latDat->GetSite(15 - delta, 15, 15).GetForce().GetMagnitudeSquared() < 1e-8);
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellCellInteractionWithGridTests);
    }
  }
}

#endif  // ONCE
