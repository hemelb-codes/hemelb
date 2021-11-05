// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <memory>

#include <catch2/catch.hpp>

#include "SimulationMaster.h"

#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/LaddFail.h"

namespace hemelb
{
  namespace tests
  {

    TEST_CASE_METHOD(helpers::FolderTestFixture, "SimulationMaster") {
      const int argc = 7;
      const char* argv[] = {
	"hemelb",
	"-in",
	"four_cube.xml",
	"-i",
	"1",
	"-ss",
	"1111"
      };

      CopyResourceToTempdir("four_cube.xml");
      CopyResourceToTempdir("four_cube.gmy");

      auto options = std::make_unique<hemelb::configuration::CommandLine>(argc, argv);
      auto master = std::make_unique<SimulationMaster<>>(*options, Comms());

      SECTION("Running a simulation creates outputs") {
	// TODO: This test is fatal if run with LADDIOLET. See ticket #605.
	LADD_FAIL();
	master->RunSimulation();
	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
	AssertPresent("results/Extracted/wholegeometryvelocityandstress.dat");
	AssertPresent("results/Extracted/centrelinepressure.dat");
	AssertPresent("results/Extracted/centrelineshearrate.dat");
	AssertPresent("results/Extracted/surfaceshearstress.dat");
	AssertPresent("results/Extracted/surfacetraction.dat");
      }
    }

  }
}
