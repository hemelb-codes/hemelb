// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_UNITTESTS_REDBLOOD_FEDOSOVVALIDATIONTESTS_H
#define HEMELB_UNITTESTS_REDBLOOD_FEDOSOVVALIDATIONTESTS_H

#include <cppunit/extensions/HelperMacros.h>
#include <boost/uuid/uuid_io.hpp>
#include <memory>

#include "SimulationMaster.h"
#include "lb/BuildSystemInterface.h"
#include "lb/lattices/D3Q19.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "redblood/stencil.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class FedosovValidationTests : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (FedosovValidationTests);
          CPPUNIT_TEST (testIntegration);CPPUNIT_TEST_SUITE_END();

          typedef Traits<>::Reinstantiate<lb::lattices::D3Q19, lb::GuoForcingLBGK>::Type Traits;
          typedef hemelb::redblood::CellController<Traits> CellControl;
          typedef SimulationMaster<Traits> MasterSim;

        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            CopyResourceToTempdir("fedosov1c.xml");
            CopyResourceToTempdir("fedosov1c.gmy");
            CopyResourceToTempdir("rbc_ico_720.msh");
            CopyResourceToTempdir("rbc_ico_1280.msh");

            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "fedosov1c.xml";
            argv[3] = "-i";
            argv[4] = "0";
            argv[5] = "-ss";
            argv[6] = "1111";
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);

            master = std::make_shared<MasterSim>(*options, Comms());
          }

          void testIntegration()
          {
            CPPUNIT_ASSERT(master);
            auto const & converter = master->GetUnitConverter();
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            CPPUNIT_ASSERT(controller);
            controller->AddCellChangeListener([&converter](const hemelb::redblood::CellContainer &cells)
            {
              static int iter = 0;
              if(cells.empty())
              {
                return;
              }
              auto cell = *cells.begin();
              if(iter % 1000 == 0)
              {
                std::stringstream filename;
                filename << cell->GetTag() << "_t_" << iter << ".vtp";
                writeVTKMesh(filename.str(), cell, converter);
              }
              ++iter;
            });

            // run the simulation
            master->RunSimulation();
            master->Finalise();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];

      };

    //CPPUNIT_TEST_SUITE_REGISTRATION (FedosovValidationTests);
    }// namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
