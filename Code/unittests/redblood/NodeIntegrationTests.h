//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLINTEGRATION_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLINTEGRATION_H

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
      //! Fake cell that contains a single node
      class NodeCell : public hemelb::redblood::CellBase
      {
        public:
          NodeCell(LatticePosition const&position)
            : hemelb::redblood::CellBase(
                {position},
                hemelb::redblood::Mesh(
                  std::make_shared<hemelb::redblood::MeshData>(
                    hemelb::redblood::MeshData{{position}, {}}
                  ),
                  std::make_shared<hemelb::redblood::MeshTopology>()
                ),
                1e0, "nope"
              )
          {
          }

          LatticeEnergy operator()() const override
          {
            return 0e0;
          }
          LatticeEnergy operator()(std::vector<LatticeForceVector> &) const override
          {
            return 0e0;
          }
          std::unique_ptr<CellBase> cloneImpl() const override
          {
            return std::unique_ptr<NodeCell>{new NodeCell(GetVertices()[0])};
          }
      };

      class NodeIntegrationTests : public helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (NodeIntegrationTests);
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::FourPoint>);
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::CosineApprox>);
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::ThreePoint>);
          CPPUNIT_TEST (testNodeWall<hemelb::redblood::stencil::TwoPoint>);
          CPPUNIT_TEST_SUITE_END();

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
            std::shared_ptr<SimulationMaster<TRAITS>> simulationMaster(
              size_t steps, Dimensionless cell, Dimensionless wall) const;

          //! Runs simulation with a single node
          template<class STENCIL>
          LatticePosition testNodeWall(Dimensionless intensity, LatticePosition where);

          //! Checks that 1-node cell moves away from the wall
          template<class STENCIL> void testNodeWall()
          {
            // no interaction because far from wall
            auto const farFromWall = testNodeWall<STENCIL>(2e0, {7e0, 7e0, 16.1e0});
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
          if(util::DoesDirectoryExist("results"))
          {
            system("rm -rf results");
          }
          DeleteXMLInput(resource, {"inlets", "inlet", "insertcell"});
          DeleteXMLInput(resource, {"inlets", "inlet", "flowextension"});
          DeleteXMLInput(resource, {"outlets", "outlet", "flowextension"});
          DeleteXMLInput(resource, {"properties", "propertyoutput"});
          ModifyXMLInput(resource, {"inlets", "inlet", "condition", "mean", "value"}, 0);
          ModifyXMLInput(resource, {"simulation", "steps", "value"}, steps);
          ModifyXMLInput(resource, {"redbloodcells", "cell2Cell", "intensity", "value"}, cell);
          ModifyXMLInput(resource, {"redbloodcells", "cell2Cell", "cutoff", "value"}, 2);
          ModifyXMLInput(resource, {"redbloodcells", "cell2Wall", "intensity", "value"}, wall);
          ModifyXMLInput(resource, {"redbloodcells", "cell2Well", "cutoff", "value"}, 2);
          auto options = std::make_shared<configuration::CommandLine>(argc, argv);
          auto const master = std::make_shared<SimulationMaster<TRAITS>>(*options, Comms());
          helpers::LatticeDataAccess(&master->GetLatticeData()).ZeroOutForces();
          return master;
        }

      template<class STENCIL>
        LatticePosition NodeIntegrationTests::testNodeWall(
            Dimensionless intensity, LatticePosition where)
        {
          typedef typename
            Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type::ChangeStencil<STENCIL>::Type
            Traits;
          typedef hemelb::redblood::CellController<Traits> CellController;

          auto const master = simulationMaster<Traits>(3, 0, intensity);
          auto controller = std::static_pointer_cast<CellController>(master->GetCellController());
          controller->AddCell(std::make_shared<NodeCell>(where));

          master->RegisterActor(*controller, 1);
          master->RunSimulation();
          master->Finalise();
          return (*controller->GetCells().begin())->GetVertices().front();
        }

      CPPUNIT_TEST_SUITE_REGISTRATION (NodeIntegrationTests);
    }
  }
}

#endif  // ONCE
