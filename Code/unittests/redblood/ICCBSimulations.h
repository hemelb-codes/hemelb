// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_UNITTESTS_REDBLOOD_ICCBSIMULATIONS_H
#define HEMELB_UNITTESTS_REDBLOOD_ICCBSIMULATIONS_H

#include <cppunit/extensions/HelperMacros.h>
#include <memory>

#include "lb/BuildSystemInterface.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "unittests/helpers/FolderTestFixture.h"
#include <boost/uuid/uuid_io.hpp>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class ICCBSimulations : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (ICCBSimulations);
          //CPPUNIT_TEST (testCoarseNetwork);
          CPPUNIT_TEST (testFineNetwork);CPPUNIT_TEST_SUITE_END();

          typedef Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type Traits;
          typedef hemelb::redblood::CellController<Traits> CellControl;
          typedef SimulationMaster<Traits> MasterSim;

        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            CopyResourceToTempdir("iccb_capillary_network.xml");
            CopyResourceToTempdir("P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_388889_1000_3.gmy");
            CopyResourceToTempdir("P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_875_1000_3.gmy");
            CopyResourceToTempdir("rbc_ico_1280.msh");

            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "iccb_capillary_network.xml";
            argv[3] = "-i";
            argv[4] = "0";
            argv[5] = "-ss";
            argv[6] = "1111";
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);
          }

          void tearDown()
          {
            master->Finalise();
            master.reset();
          }

          void testCoarseNetwork()
          {
            ModifyXMLInput("iccb_capillary_network.xml",
                           { "simulation", "step_length", "value" },
                           1.488715e-07);
            ModifyXMLInput("iccb_capillary_network.xml",
                           { "simulation", "voxel_size", "value" },
                           8.75e-07);
            ModifyXMLInput("iccb_capillary_network.xml",
                           { "simulation", "origin", "value" },
                           "(1.2574028194e-05,9.44294101e-06,-1.13739131093e-05)");
            ModifyXMLInput("iccb_capillary_network.xml",
                           { "geometry", "datafile", "path" },
                           "P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_875_1000_3.gmy");

            ModifyXMLInput("iccb_capillary_network.xml", { "redbloodcells",
                                                           "cells",
                                                           "cell",
                                                           "moduli",
                                                           "bending",
                                                           "value" },
                           1.135e-05);
            ModifyXMLInput("iccb_capillary_network.xml", { "redbloodcells",
                                                           "cells",
                                                           "cell",
                                                           "moduli",
                                                           "strain",
                                                           "value" },
                           0.0004597);
            ModifyXMLInput("iccb_capillary_network.xml", { "redbloodcells",
                                                           "cell2Cell",
                                                           "intensity",
                                                           "value" },
                           1.135e-05);
            ModifyXMLInput("iccb_capillary_network.xml", { "redbloodcells",
                                                           "cell2Wall",
                                                           "intensity",
                                                           "value" },
                           1.135e-05);

            // The MasterSim object has to be created after the XML file has been modified for those changes to be taken into account
            master = std::make_shared<MasterSim>(*options, Comms());
            CPPUNIT_ASSERT(master);
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            CPPUNIT_ASSERT(controller);

            unsigned timestep = 0;
            auto output_callback =
                [this, &timestep](const hemelb::redblood::CellContainer & cells)
                {
                  if ((timestep % 1000) == 0)
                  {
                    for (auto cell: cells)
                    {
                      std::stringstream filename;
                      filename << cell->GetTag() << "_t_" << timestep << ".vtp";
                      hemelb::redblood::writeVTKMesh(filename.str(), cell, this->master->GetUnitConverter());
                    }
                  }
                  timestep++;
                };
            controller->AddCellChangeListener(output_callback);

            // run the simulation
            master->RunSimulation();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

          void testFineNetwork()
          {
            master = std::make_shared<MasterSim>(*options, Comms());
            CPPUNIT_ASSERT(master);
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            CPPUNIT_ASSERT(controller);

            unsigned timestep = 0;
            auto output_callback =
                [this, &timestep](const hemelb::redblood::CellContainer & cells)
                {
                  if ((timestep % 1000) == 0)
                  {
                    for (auto cell: cells)
                    {
                      std::stringstream filename;
                      filename << cell->GetTag() << "_t_" << timestep << ".vtp";
                      hemelb::redblood::writeVTKMesh(filename.str(), cell, this->master->GetUnitConverter());
                    }
                  }
                  timestep++;
                };
            controller->AddCellChangeListener(output_callback);

            // run the simulation
            master->RunSimulation();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];
      };

    // Don't register the unit test so it's not run by default as part of CI.
    // Uncomment the line below in order to run the test with: ./unittests_hemelb hemelb::unittests::redblood::ICCBSimulations::testSimpleTube
    //CPPUNIT_TEST_SUITE_REGISTRATION (ICCBSimulations);

    }// namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
