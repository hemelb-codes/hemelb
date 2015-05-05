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
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/helpers/LatticeDataAccess.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellIntegrationTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (CellIntegrationTests);
            CPPUNIT_TEST (testCellOutOfBounds);
            CPPUNIT_TEST (testIntegration);
            CPPUNIT_TEST (testIntegrationWithoutCells);
            CPPUNIT_TEST (testSwamped);
          CPPUNIT_TEST_SUITE_END();

          typedef Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type Traits;
          typedef CellController<Traits::Kernel> CellControll;
          typedef SimulationMaster<Traits> MasterSim;
        public:
          void setUp()
          {
            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "large_cylinder.xml";
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-ss";
            argv[6] = "1111";
            FolderTestFixture::setUp();
            CopyResourceToTempdir("red_blood_cell.txt");
            TiXmlDocument doc(resources::Resource("large_cylinder.xml").Path());
            CopyResourceToTempdir("large_cylinder.xml");
            std::vector<std::string> intel;
            intel.push_back("simulation"); intel.push_back("steps"); intel.push_back("value");
            ModifyXMLInput("large_cylinder.xml", std::move(intel), 8);
            CopyResourceToTempdir("large_cylinder.gmy");
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);
 
            auto cell = std::make_shared<Cell>(icoSphere(4));
            cell->moduli = Cell::Moduli(1e-6, 1e-6, 1e-6, 1e-6);
            cells.insert(cell);
 
            master = std::make_shared<MasterSim>(*options, Comms());
            helpers::LatticeDataAccess(&master->GetLatticeData()).ZeroOutForces();
          }

          void tearDown()
          {
            master->Finalise();
            master.reset();
          }

          // No errors when interpolation/spreading hits nodes out of bounds
          void testCellOutOfBounds()
          {
            (*cells.begin())->operator+=(master->GetLatticeData().GetGlobalSiteMins() * 2);
            auto controller = std::make_shared<CellControll>(master->GetLatticeData(), cells);
            master->RegisterActor(*controller, 1);
            master->RunSimulation();
            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }
  
          // Check that the particles move and result in some force acting on the fluid
          void testIntegration()
          {
            // setup cell position
            auto const &latticeData = master->GetLatticeData();
            auto const mid = LatticePosition(latticeData.GetGlobalSiteMaxes()
                  + latticeData.GetGlobalSiteMins()) * 0.5;
            (**cells.begin()) += mid - (*cells.begin())->GetBarycenter();
            (**cells.begin()) += LatticePosition(0, 0, 8 - mid.z);
            (**cells.begin()) *= 5.0;
            auto controller = std::make_shared<CellControll>(master->GetLatticeData(), cells);
            auto const barycenter = (*cells.begin())->GetBarycenter();

            // run
            master->RegisterActor(*controller, 1);
            master->RunSimulation();

            // check position of cell has changed
            auto const moved = (*cells.begin())->GetBarycenter();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.x - moved.x, 0e0, 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.y - moved.y, 0e0, 1e-6);
            CPPUNIT_ASSERT(std::abs(barycenter.z - moved.z) > 1e-3);

            // check there is force on one of the lattice site near a node
            // node position is guessed at from geometry
            auto const nodepos = mid + LatticePosition(0, 0, 8 - 5 - mid.z);
            auto const force = latticeData.GetSite(nodepos).GetForce() ;
            CPPUNIT_ASSERT(std::abs(force.z) > 1e-4);

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

          // Check that the particles move and result in some force acting on the fluid
          void testIntegrationWithoutCells()
          {
            // setup cell position
            cells.clear();
            auto controller = std::make_shared<CellControll>(master->GetLatticeData(), cells);

            // run
            master->RegisterActor(*controller, 1);
            master->RunSimulation();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

          // If particle is rigid, then  it should not move within the length of the calculation
          // This test is meaningfull only if the other two are not failing
          void testSwamped()
          {
            // setup cell position
            auto const &latticeData = master->GetLatticeData();
            auto const mid = LatticePosition(latticeData.GetGlobalSiteMaxes()
                  + latticeData.GetGlobalSiteMins()) * 0.5;
            (*(*cells.begin())) += mid - (*cells.begin())->GetBarycenter();
            (*(*cells.begin())) += LatticePosition(0, 0, 8 - mid.z);
            (*(*cells.begin())) *= 5.0;
            auto controller = std::make_shared<CellControll>(master->GetLatticeData(), cells);
            auto const barycenter = (*cells.begin())->GetBarycenter();
            std::dynamic_pointer_cast<Cell>(*cells.begin())->moduli
                = Cell::Moduli(1e0, 1e0, 1e0, 1e0, 1e0);

            // run
            master->RegisterActor(*controller, 1);
            master->RunSimulation();

            // check position of cell has *not* changed
            auto const moved = (*cells.begin())->GetBarycenter();
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.x - moved.x, 0e0, 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.y - moved.y, 0e0, 1e-6);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(barycenter.z - moved.z, 0e0, 1e-6);

            // check that force on lattice is very large
            auto const nodepos = mid + LatticePosition(0, 0, 8 - 5 - mid.z);
            auto const force = latticeData.GetSite(nodepos).GetForce() ;
            CPPUNIT_ASSERT(std::abs(force.z) > 1e2);

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          hemelb::redblood::CellContainer cells;
          int const argc = 7;
          char const * argv [7];
      };


      CPPUNIT_TEST_SUITE_REGISTRATION (CellIntegrationTests);
    }
  }
}

#endif  // ONCE
