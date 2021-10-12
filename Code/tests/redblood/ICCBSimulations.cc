// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>
#include <boost/uuid/uuid_io.hpp>
#include <catch2/catch.hpp>

#include "lb/BuildSystemInterface.h"
#include "Traits.h"
#include "SimulationMaster.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    class ICCBSimulations : public hemelb::tests::helpers::FolderTestFixture
    {
      typedef Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type Traits;
      typedef hemelb::redblood::CellController<Traits> CellControl;
      typedef SimulationMaster<Traits> MasterSim;
    public:
      ICCBSimulations() : FolderTestFixture() {
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

      ~ICCBSimulations() {
	master->Finalise();
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

	// The MasterSim object has to be created after the XML file
	// has been modified for those changes to be taken into
	// account
	master = std::make_shared<MasterSim>(*options, Comms());
	REQUIRE(master);
	auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
	REQUIRE(controller);

	unsigned timestep = 0;
	auto output_callback =
	  [this, &timestep](const hemelb::redblood::CellContainer & cells) {
	  if ((timestep % 1000) == 0) {
	    for (auto& cell: cells) {
	      std::stringstream filename;
	      filename << cell->GetTag() << "_t_" << timestep << ".vtp";
	      io.writeFile(filename.str(), *cell, this->master->GetUnitConverter());
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

      void testFineNetwork() {
	master = std::make_shared<MasterSim>(*options, Comms());
	REQUIRE(master);
	auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
	REQUIRE(controller);

	unsigned timestep = 0;
	auto output_callback =
	  [this, &timestep](const hemelb::redblood::CellContainer & cells) {
	  if ((timestep % 1000) == 0) {
	    for (auto& cell: cells) {
	      std::stringstream filename;
	      filename << cell->GetTag() << "_t_" << timestep << ".vtp";
	      io.writeFile(filename.str(), *cell, this->master->GetUnitConverter());
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
      redblood::VTKMeshIO io = {};
    };

    // These tests were never run in CI and are missing inputs.
    METHOD_AS_TEST_CASE(ICCBSimulations::testCoarseNetwork,
			"Some simulation on a coarser mesh",
			"[redblood][.broken]");
    METHOD_AS_TEST_CASE(ICCBSimulations::testFineNetwork,
			"Some simulation on a finer mesh",
			"[redblood][.broken]");

  } // namespace tests
} // namespace hemelb

