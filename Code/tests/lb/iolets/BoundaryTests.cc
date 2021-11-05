// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "resources/Resource.h"
#include "lb/iolets/BoundaryValues.h"
#include "tests/helpers/LaddFail.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::lb::iolets;
    using namespace resources;

    // Class asserting behaviour of boundary-collection objects and
    // the boundaries (in- and out- lets) within them.
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "BoundaryTests") {
      auto pressureToDensity = [&](double pressure) -> double {
	double inverseVelocity = simConfig->GetTimeStepLength() / simConfig->GetVoxelSize();
	return 1
	  + pressure * mmHg_TO_PASCAL * inverseVelocity * inverseVelocity
	  / (Cs2 * lbmParams->GetFluidDensity());
      };

      auto inlets = std::make_unique<BoundaryValues>(
	  geometry::INLET_TYPE,
	  latDat,
	  simConfig->GetInlets(),
	  simState.get(),
	  Comms(),
	  *unitConverter);

      SECTION("TestConstruct") {
	double targetStartDensity = unitConverter->ConvertPressureToLatticeUnits(80.0 - 1.0) / Cs2;

	REQUIRE(Approx(targetStartDensity) == inlets->GetBoundaryDensity(0));
      }

      SECTION("TestUpdate") {
	REQUIRE(Approx(pressureToDensity(80.0 - 1.0)) == inlets->GetBoundaryDensity(0));

	while (simState->Get0IndexedTimeStep() < simState->GetTotalTimeSteps() / 20) {
	  simState->Increment();
	}

	REQUIRE(Approx(pressureToDensity(80.0 + 1.0)) == inlets->GetBoundaryDensity(0));
	
	while (simState->Get0IndexedTimeStep() < simState->GetTotalTimeSteps() / 10) {
	  simState->Increment();
	}

	REQUIRE(Approx(pressureToDensity(80.0 - 1.0)) == inlets->GetBoundaryDensity(0));
      }

      SECTION("TestUpdateFile") {
	LADD_FAIL();
	CopyResourceToTempdir("iolet.txt");
	MoveToTempdir();

	auto fileInletConfig = std::unique_ptr<configuration::SimConfig>{configuration::SimConfig::New(Resource("config_file_inlet.xml").Path())};

	// Reloading simState to ensure the time step size from config_file_inlet.xml is indeed used in HemeLB.
	simState.reset(new hemelb::lb::SimulationState{fileInletConfig->GetTimeStepLength(), fileInletConfig->GetTotalTimeSteps()});

	inlets.reset(new BoundaryValues(hemelb::geometry::INLET_TYPE,
					latDat,
					fileInletConfig->GetInlets(),
					simState.get(),
					Comms(),
					*unitConverter));

	REQUIRE(Approx(pressureToDensity(78.0)) == inlets->GetBoundaryDensity(0));

	while (simState->Get0IndexedTimeStep() < simState->GetTotalTimeSteps() / 2) {
	  simState->Increment();
	}

	REQUIRE(Approx(pressureToDensity(82.0)) == inlets->GetBoundaryDensity(0));

	while (simState->Get0IndexedTimeStep() < simState->GetTotalTimeSteps()) {
	  simState->Increment();
	}

	REQUIRE(Approx(pressureToDensity(78.0)) == inlets->GetBoundaryDensity(0));
      }
    }
  }
}
