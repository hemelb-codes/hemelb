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

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class CellIntegrationTests : public FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (CellIntegrationTests);
          CPPUNIT_TEST (testIntegration);
          CPPUNIT_TEST_SUITE_END();

          typedef lb::lattices::HEMELB_LATTICE Lattice;
          typedef lb::HEMELB_KERNEL<Lattice>::Type Kernel;
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
            CopyResourceToTempdir("large_cylinder.xml");
            CopyResourceToTempdir("large_cylinder.gmy");
            try
            {
              options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);
            }
            catch (hemelb::io::xml::ChildError& e)
            {
              std::cout << e.what() << std::endl;
              throw;
            }

            auto cell = std::make_shared<Cell>(icoSphere(4));
            *cell += hemelb::LatticePosition(25, 25, 25);
            cells.emplace_back(cell);
          }

          void testIntegration()
          {
            auto master = std::make_shared<SimulationMaster>(*options, Comms());
            auto controller
              = std::make_shared<CellController<Kernel>>(master->GetLatticeData(), cells);
            master->SetCellController(controller);
            master->RunSimulation();
            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }
        private:
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
