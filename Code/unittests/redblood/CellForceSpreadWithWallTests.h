//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELL_FORCE_SPREAD_WITH_WALL_TESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_CELL_FORCE_SPREAD_WITH_WALL_TESTS_H

#include <cppunit/TestFixture.h>
#include "redblood/GridAndCell.h"
#include "redblood/Facet.impl.cc"
#include "unittests/redblood/Fixtures.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellForceSpreadWithWallTests : public SquareDuctTetrahedronFixture
      {
          CPPUNIT_TEST_SUITE (CellForceSpreadWithWallTests);
          CPPUNIT_TEST (testAPIAssumptions);
          CPPUNIT_TEST (testNode2WallCutoff<stencil::FourPoint>);
          CPPUNIT_TEST (testNode2WallCutoff<stencil::CosineApprox>);
          CPPUNIT_TEST (testNode2WallCutoff<stencil::ThreePoint>);
          CPPUNIT_TEST (testNode2WallCutoff<stencil::TwoPoint>);
          CPPUNIT_TEST_SUITE_END();
          typedef lb::lattices::D3Q15 D3Q15;
          typedef lb::kernels::LBGK<D3Q15> Kernel;

        public:
          void testAPIAssumptions();
          template<class STENCIL> void testNode2WallCutoff();

        protected:
          int siteID(LatticeVector const &) const;
          LatticeVector GetSolidWall() const;
          virtual Mesh initial_mesh() const
          {
            return triangleMesh();
          }
          virtual size_t refinement() const
          {
            return 0;
          }
      };

      // IsValidLatticeSite doesn't always return false when site is invalid
      // It seems to be more of a quick and innacurate check
      int CellForceSpreadWithWallTests::siteID(LatticeVector const &position) const
      {
        site_t siteid;
        proc_t procid;
        bool const valid = latDat->GetContiguousSiteId(position, procid, siteid);
        return valid ?
          siteid :
          -1;
      }
      LatticeVector CellForceSpreadWithWallTests::GetSolidWall() const
      {
        return latDat->GetGlobalSiteMaxes() - LatticeVector(CubeSize() / 2, -1, CubeSize() / 2);
      }

      void CellForceSpreadWithWallTests::testAPIAssumptions()
      {
        // Nearest solid node. Does not have associated data.
        LatticeVector const solid = GetSolidWall();
        // Nearest wet wall node. Has fluid particles and is closest node to wall.
        LatticeVector const wetwall = solid + LatticeVector(0, -1, 0);
        LatticeVector const left = wetwall + LatticeVector(-1, 0, 0);
        LatticeVector const right = wetwall + LatticeVector(1, 0, 0);
        LatticeVector const fluid = wetwall + LatticeVector(0, -1, 0);

        // Checking I know where walls and fluid sites are.
        // This is used subsequently in the tests
        CPPUNIT_ASSERT(siteID(solid) < 0);
        CPPUNIT_ASSERT(siteID(wetwall) >= 0);
        CPPUNIT_ASSERT(siteID(left) >= 0);
        CPPUNIT_ASSERT(siteID(right) >= 0);
        CPPUNIT_ASSERT(not latDat->GetSite(siteID(fluid)).IsSolid());
        CPPUNIT_ASSERT(not latDat->GetSite(siteID(wetwall)).IsSolid());
        CPPUNIT_ASSERT(not latDat->GetSite(siteID(left)).IsSolid());
        CPPUNIT_ASSERT(not latDat->GetSite(siteID(right)).IsSolid());
        CPPUNIT_ASSERT(not latDat->GetSite(siteID(fluid)).IsWall());
        CPPUNIT_ASSERT(latDat->GetSite(siteID(wetwall)).IsWall());
        CPPUNIT_ASSERT(latDat->GetSite(siteID(left)).IsWall());
        CPPUNIT_ASSERT(latDat->GetSite(siteID(right)).IsWall());
      }

      template<class STENCIL> void CellForceSpreadWithWallTests::testNode2WallCutoff()
      {
        // Fluid sites next to wall
        LatticeVector const solid = GetSolidWall();
        LatticeVector const wetwall = solid + LatticeVector(0, -1, 0);
        LatticeVector const left = wetwall + LatticeVector(-1, 0, 0);
        LatticeVector const right = wetwall + LatticeVector(1, 0, 0);

        // Makes sure there are no force except node-Wall interaction
        helpers::ZeroOutForces(latDat);
        mesh.moduli = Cell::Moduli();
        mesh.nodeWall.intensity = 1.;
        mesh += LatticePosition(10, 10, 10);
        // Min distance to wall, so we can check cutoff
        helpers::SetWallDistance(latDat, 0.3);

        struct TestCase
        {
           LatticeDistance cutoff;
           LatticePosition position;
           bool atWall;
           bool atLeft;
           bool atRight;
        };
        std::vector<TestCase> const testcases
        {
          {0.25, wetwall.cast<LatticeDistance>(), false, false, false},
          {0.35, wetwall.cast<LatticeDistance>(), true, false, false},
          {0.35, wetwall.cast<LatticeDistance>() - LatticePosition(0, 0.2, 0), false, false, false},
          {
            std::sqrt(0.2 * 0.2 + 1e0*1e0) + 0.1,
            wetwall.cast<LatticeDistance>() - LatticePosition(0, 0.2, 0),
            true,
            Dimensionless(STENCIL::GetRange()) * 0.5 > std::sqrt(0.2 * 0.2 + 1e0*1e0) + 0.1,
            Dimensionless(STENCIL::GetRange()) * 0.5 > std::sqrt(0.2 * 0.2 + 1e0*1e0) + 0.1
          },
        };

        // Loop over test cases breaks on special marker (negative) cutoff
        for (auto const &testcase: testcases)
        {
          helpers::ZeroOutForces(latDat);
          mesh.GetVertices()[0] = testcase.position;
          mesh.nodeWall.cutoff = testcase.cutoff;

          // Shared pointer does not deallocate. It's a fake so we still use forcesOnGrid interface.
          std::shared_ptr<CellBase> ptr_mesh(&mesh, [](CellBase*){});
          forcesOnGrid < D3Q15, STENCIL > (ptr_mesh, *latDat);

          bool const atWall = not helpers::is_zero(latDat->GetSite(wetwall).GetForce());
          bool const atLeft = not helpers::is_zero(latDat->GetSite(left).GetForce());
          bool const atRight = not helpers::is_zero(latDat->GetSite(right).GetForce());
          CPPUNIT_ASSERT_EQUAL(testcase.atWall, atWall);
          CPPUNIT_ASSERT_EQUAL(testcase.atRight, atRight);
          CPPUNIT_ASSERT_EQUAL(testcase.atLeft, atLeft);
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellForceSpreadWithWallTests);
    }
  }
}

#endif  // ONCE
