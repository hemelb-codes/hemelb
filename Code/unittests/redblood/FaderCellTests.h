//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_FADERCELLS_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_FADERCELLS_TESTS_H

#include <cppunit/TestFixture.h>

#include "unittests/redblood/Fixtures.h"
#include "redblood/FaderCell.h"
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
        CPPUNIT_TEST_SUITE(FaderCellTests);
          CPPUNIT_TEST(testFade);
        CPPUNIT_TEST_SUITE_END();
        public:
          void testFade();
      };

      class FaderCellTests::DummyCell : public CellBase
      {
        public:
          DummyCell() : CellBase(icoSphere())
          {
          }

          PhysicalEnergy operator()() const override
          {
            return 1.0;
          }
          PhysicalEnergy operator()(std::vector<LatticeForceVector> & forces) const override
          {
            int i(1);
            for(auto& force: forces)
            {
              force = LatticeForceVector(0, i, 8);
              ++i;
            }
            return 1e0;
          }
          LatticeForceVector WallInteractionForce(
              LatticePosition const& vertex, LatticePosition const &wall) const override
          {
            return {1, 2, 3};
          }
      };

      void FaderCellTests::testFade()
      {
        // normals of inlet and outlet face in different direction
        FlowExtension const inlet (util::Vector3D<Dimensionless>(-1, 0, 0), LatticePosition(3, 2, 2), 2.0, 2, 1.8);
        FlowExtension const outlet(util::Vector3D<Dimensionless>( 1, 0, 0), LatticePosition(8, 2, 2), 2.0, 2, 1.8);
        auto dummyCell = std::make_shared<FaderCellTests::DummyCell>();
        std::vector<FlowExtension> extensions; extensions.push_back(inlet); extensions.push_back(outlet);
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
        auto force = dummyCell->WallInteractionForce(zero, zero);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, force.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2e0, force.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3e0, force.z, 1e-8);

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
        // wall interactions do *not* fade
        force = dummyCell->WallInteractionForce(zero, zero);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, force.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2e0, force.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3e0, force.z, 1e-8);

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
        force = dummyCell->WallInteractionForce(zero, zero);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1e0, force.x, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2e0, force.y, 1e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3e0, force.z, 1e-8);
      }




      CPPUNIT_TEST_SUITE_REGISTRATION (FaderCellTests);
    }
  }
} // hemelb::redblood
#endif
