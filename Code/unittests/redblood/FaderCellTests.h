// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_FADERCELLTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_FADERCELLTESTS_H

#include <cppunit/TestFixture.h>

#include "unittests/redblood/Fixtures.h"
#include "redblood/FaderCell.h"
#include "redblood/CellCell.h"
#include "redblood/FlowExtension.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class FaderCellTests : public CppUnit::TestFixture
      {
          class DummyCell;
          CPPUNIT_TEST_SUITE (FaderCellTests);
          CPPUNIT_TEST (testFade);
          CPPUNIT_TEST (testDivideAndConquer);CPPUNIT_TEST_SUITE_END();
        public:
          // Cell does fade
          void testFade();
          // Divide and conquer still finds neighboring node of a fader cell
          void testDivideAndConquer();
      };

      class FaderCellTests::DummyCell : public CellBase
      {
        public:
          DummyCell() :
              CellBase(icoSphere())
          {
          }

          LatticeEnergy operator()() const override
          {
            return 1.0;
          }
          LatticeEnergy operator()(std::vector<LatticeForceVector> & forces) const override
          {
            int i(1);
            for (auto& force : forces)
            {
              force = LatticeForceVector(0, i, 8);
              ++i;
            }
            return 1e0;
          }
        private:
          std::unique_ptr<CellBase> cloneImpl() const override
          {
            return std::unique_ptr<FaderCellTests::DummyCell>(new DummyCell());
          }
      };

      void FaderCellTests::testFade()
      {
        // normals of inlet and outlet face in different direction
        FlowExtension const inlet(util::Vector3D<Dimensionless>(-1, 0, 0),
                                  LatticePosition(3, 2, 2),
                                  2.0,
                                  2,
                                  1.8);
        FlowExtension const outlet(util::Vector3D<Dimensionless>(1, 0, 0),
                                   LatticePosition(8, 2, 2),
                                   2.0,
                                   2,
                                   1.8);
        auto dummyCell = std::make_shared<FaderCellTests::DummyCell>();
        std::vector<FlowExtension> extensions;
        extensions.push_back(inlet);
        extensions.push_back(outlet);
        auto const fadingCell = std::make_shared<FaderCell>(dummyCell, extensions);

        const LatticePosition zero(0, 0, 0);
        std::vector<LatticeForceVector> forces(2, zero);
        // No fading here
        *fadingCell += LatticePosition(5, 1, 1) - fadingCell->GetBarycenter();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, fadingCell->Energy(forces), 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, forces.front().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, forces.front().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(8e0, forces.front().z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0e0, forces.back().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2e0, forces.back().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(8e0, forces.back().z, 1e-8);

        // Now fades to 0.7
        forces = std::vector<LatticeForceVector>(2, zero);
        *fadingCell += LatticePosition(3.0 - inlet.fadeLength * 0.3, 1, 1)
            - fadingCell->GetBarycenter();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 1e0, fadingCell->Energy(forces), 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 0e0, forces.front().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 1e0, forces.front().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 8e0, forces.front().z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 0e0, forces.back().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 2e0, forces.back().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 8e0, forces.back().z, 1e-8);

        // fades to 0.7 in outlet
        forces = std::vector<LatticeForceVector>(2, zero);
        *fadingCell += LatticePosition(8.0 + inlet.fadeLength * 0.3, 1, 1)
            - fadingCell->GetBarycenter();
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 1e0, fadingCell->Energy(forces), 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 0e0, forces.front().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 1e0, forces.front().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 8e0, forces.front().z, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 0e0, forces.back().x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 2e0, forces.back().y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7 * 8e0, forces.back().z, 1e-8);
      }

      void FaderCellTests::testDivideAndConquer()
      {
        FlowExtension const inlet(util::Vector3D<Dimensionless>(-1, 0, 0),
                                  LatticePosition(3, 2, 2),
                                  2.0,
                                  2,
                                  1.8);
        FlowExtension const outlet(util::Vector3D<Dimensionless>(1, 0, 0),
                                   LatticePosition(8, 2, 2),
                                   2.0,
                                   2,
                                   1.8);

        auto const cell = std::make_shared<Cell>(tetrahedron());
        *cell *= 5e0;
        auto const fader0 =
            std::make_shared<FaderCell>(cell, std::vector<FlowExtension> { inlet, outlet });

        // Check cloning while we are at it
        std::shared_ptr<FaderCell> const fader1(fader0->clone().release());
        CPPUNIT_ASSERT(fader0 != fader1);
        CPPUNIT_ASSERT(&fader0->GetVertices() != &fader1->GetVertices());
        CPPUNIT_ASSERT(fader0->GetTag() != fader1->GetTag());

        // Now move second cell to have a close node to another
        auto const trans = fader0->GetVertices().front() - fader0->GetVertices().back();
        *fader1 += trans + trans.GetNormalised() * 0.5;

        // Add them to a divide and conquer object
        DivideConquerCells dnc( { fader0, fader1 }, 100e0, 1e0);

        auto range = dnc.pair_begin(0.6);
        CPPUNIT_ASSERT(range.is_valid());
        CPPUNIT_ASSERT(not (++range));
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (FaderCellTests);
    }
  }
} // hemelb::redblood
#endif
