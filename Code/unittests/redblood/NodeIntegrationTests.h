// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REDBLOOD_NODEINTEGRATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_NODEINTEGRATIONTESTS_H

#include <cppunit/TestFixture.h>
#include <cstdio>
#include <unistd.h>
#include "util/fileutils.h"
#include "Traits.h"
#include "SimulationMaster.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class NodeIntegrationTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (NodeIntegrationTests);
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::FourPoint> );
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::CosineApprox> );
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::ThreePoint> );
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::TwoPoint> );
          CPPUNIT_TEST (testNodeNode<hemelb::redblood::stencil::FourPoint> );
          CPPUNIT_TEST (testNodeNode<hemelb::redblood::stencil::CosineApprox> );
          CPPUNIT_TEST (testNodeNode<hemelb::redblood::stencil::ThreePoint> );
          CPPUNIT_TEST (testNodeNode<hemelb::redblood::stencil::TwoPoint> );CPPUNIT_TEST_SUITE_END();

        public:
          void setUp()
          {
            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = resource.c_str();
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-ss";
            argv[6] = "1111";
            FolderTestFixture::setUp();
            CopyResourceToTempdir("red_blood_cell.txt");
            CopyResourceToTempdir("large_cylinder.gmy");
          }

          //! Creates a simulation. Does not run it.
          template<class TRAITS>
          std::shared_ptr<SimulationMaster<TRAITS>> simulationMaster(size_t steps,
                                                                     Dimensionless cell,
                                                                     Dimensionless wall) const;

          //! Runs simulation with a single node
          template<class STENCIL>
          LatticePosition testNodeWall(Dimensionless intensity, LatticePosition const &where);

          //! Checks that 1-node cell moves away from the wall
          template<class STENCIL> void testNodeWall()
          {
            // no interaction because far from wall
            auto const farFromWall = testNodeWall<STENCIL>(2e0, { 7e0, 7e0, 16.1e0 });
            CPPUNIT_ASSERT_DOUBLES_EQUAL(7e0, farFromWall.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(7e0, farFromWall.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(16.1e0, farFromWall.z, 1e-8);

            // no interaction because intensity is zero
            LatticePosition const nearWall(3.5, 2.5, 16.0);
            auto const noForce = testNodeWall<STENCIL>(0e0, nearWall);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(nearWall.x, noForce.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(nearWall.y, noForce.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(nearWall.z, noForce.z, 1e-8);

            // yes we can
            auto const withForce = testNodeWall<STENCIL>(1e0, nearWall);
            CPPUNIT_ASSERT(nearWall.x + 1e-4 < withForce.x);
            CPPUNIT_ASSERT(nearWall.y + 1e-4 < withForce.y);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(nearWall.z, withForce.z, 1e-8);

            // Forces depend on node location, so may not be symmetric
            auto const nonSym = testNodeWall<STENCIL>(1e0, nearWall + LatticePosition(0, 0, 0.1));
            CPPUNIT_ASSERT(nearWall.x + 1e-4 < nonSym.x);
            CPPUNIT_ASSERT(nearWall.y + 1e-4 < nonSym.y);
            CPPUNIT_ASSERT(nearWall.z + 1e-4 < nonSym.z);
          }

          //! Runs simulation with two nodes
          template<class STENCIL>
          std::pair<LatticePosition, LatticePosition> testNodeNode(Dimensionless intensity,
                                                                   LatticePosition const &n0,
                                                                   LatticePosition const & n1);

          //! Checks that 2 nodes move away from one another
          template<class STENCIL> void testNodeNode()
          {
            // no interaction because far from wall
            LatticePosition const center(7, 7, 20);
            LatticePosition const z(0, 0, 1);
            auto const distant = testNodeNode<STENCIL>(2e0, center - z * 2e0, center + z * 2e0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, distant.first.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, distant.first.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.z - z.z * 2, distant.first.z, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, distant.second.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, distant.second.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.z + z.z * 2, distant.second.z, 1e-8);

            // no interaction because intensity is zero
            auto const noForce = testNodeNode<STENCIL>(0e0, center - z * 0.2, center + z * 0.2);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, noForce.first.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, noForce.first.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.z - z.z * 0.2, noForce.first.z, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, noForce.second.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, noForce.second.y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.z + z.z * 0.2, noForce.second.z, 1e-8);

            // yes we can
            auto const withForce = testNodeNode<STENCIL>(1e0, center - z * 0.2, center + z * 0.2);
            HEMELB_CAPTURE2(withForce.first, withForce.second);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, withForce.first.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, withForce.first.y, 1e-8);
            CPPUNIT_ASSERT(center.z - z.z * 0.2 - 1e-4 > withForce.first.z);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.x, withForce.second.x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.y, withForce.second.y, 1e-8);
            CPPUNIT_ASSERT(center.z + z.z * 0.2 + 1e-4 < withForce.second.z);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(center.z - withForce.second.z,
                                         withForce.first.z - center.z,
                                         1e-8);
          }

        private:
          int const argc = 7;
          char const * argv[7];
          std::string const resource = "large_cylinder_rbc.xml";
      };

      template<class TRAITS>
      std::shared_ptr<SimulationMaster<TRAITS>> NodeIntegrationTests::simulationMaster(
          size_t steps, Dimensionless cell, Dimensionless wall) const
      {
        CopyResourceToTempdir(resource);
        if (util::DoesDirectoryExist("results"))
        {
          system("rm -rf results");
        }
        DeleteXMLInput(resource, { "inlets", "inlet", "insertcell" });
        DeleteXMLInput(resource, { "inlets", "inlet", "flowextension" });
        DeleteXMLInput(resource, { "outlets", "outlet", "flowextension" });
        DeleteXMLInput(resource, { "properties", "propertyoutput" });
        ModifyXMLInput(resource, { "inlets", "inlet", "condition", "mean", "value" }, 0);
        ModifyXMLInput(resource, { "simulation", "steps", "value" }, steps);
        ModifyXMLInput(resource, { "redbloodcells", "cell2Cell", "intensity", "value" }, cell);
        ModifyXMLInput(resource, { "redbloodcells", "cell2Cell", "cutoff", "value" }, 2);
        ModifyXMLInput(resource, { "redbloodcells", "cell2Wall", "intensity", "value" }, wall);
        ModifyXMLInput(resource, { "redbloodcells", "cell2Wall", "cutoff", "value" }, 2);
        auto options = std::make_shared<configuration::CommandLine>(argc, argv);
        auto const master = std::make_shared<SimulationMaster<TRAITS>>(*options, Comms());
        helpers::LatticeDataAccess(&master->GetLatticeData()).ZeroOutForces();
        return master;
      }

      template<class STENCIL>
      LatticePosition NodeIntegrationTests::testNodeWall(Dimensionless intensity,
                                                         LatticePosition const & where)
      {
        typedef typename Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type::ChangeStencil<STENCIL>::Type Traits;
        typedef hemelb::redblood::CellController<Traits> CellController;

        auto const master = simulationMaster<Traits>(3, 0, intensity);
        auto controller = std::static_pointer_cast<CellController>(master->GetCellController());
        assert(controller); // XML file contains RBC problem definition, therefore a CellController should exist already
        controller->AddCell(std::make_shared<NodeCell>(where));

        master->RunSimulation();
        return (*controller->GetCells().begin())->GetVertices().front();
      }

      template<class STENCIL>
      std::pair<LatticePosition, LatticePosition> NodeIntegrationTests::testNodeNode(
          Dimensionless intensity, LatticePosition const & n0, LatticePosition const & n1)
      {
        typedef typename Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type::ChangeStencil<STENCIL>::Type Traits;
        typedef hemelb::redblood::CellController<Traits> CellController;

        auto const master = simulationMaster<Traits>(3, intensity, 0);
        auto controller = std::static_pointer_cast<CellController>(master->GetCellController());
        assert(controller); // XML file contains RBC problem definition, therefore a CellController should exist already
        auto const firstCell = std::make_shared<NodeCell>(n0);
        auto const secondCell = std::make_shared<NodeCell>(n1);
        controller->AddCell(firstCell);
        controller->AddCell(secondCell);

        master->RunSimulation();
        return
        { firstCell->GetVertices().front(), secondCell->GetVertices().front()};
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (NodeIntegrationTests);
    }
  }
}

#endif  // ONCE
