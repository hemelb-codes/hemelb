// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iomanip>

#include <catch2/catch.hpp>
#include <ctemplate/template.h>

#include "lb/IncompressibilityChecker.h"
#include "lb/IncompressibilityChecker.hpp"
#include "reporting/BuildInfo.h"
#include "reporting/Reporter.h"
#include "reporting/Timers.h"
#include "reporting/Timers.hpp"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/HasCommsTestFixture.h"
#include "tests/lb/BroadcastMocks.h"
#include "tests/lb/LbTestsHelper.h"
#include "tests/reporting/Mocks.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::reporting;

    using TimersMock = TimersBase<ClockMock, MPICommsMock>;
    using IncompressibilityCheckerMock = lb::IncompressibilityChecker<net::BroadcastMockRootNode>;

    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "ReporterTests") {
      Reporter*reporter;

      // We need two sets of timers because the incompressibility
      // checker is not templated over timing policy.
      TimersMock *mockTimers;
      reporting::Timers* realTimers;
	  
      lb::SimulationState *state;
      lb::MacroscopicPropertyCache* cache;
      IncompressibilityCheckerMock *incompChecker;

      net::Net *net;
      hemelb::tests::FourCubeLatticeData* latticeData;
      reporting::BuildInfo *buildInfo;

      mockTimers = new TimersMock(Comms());
      realTimers = new reporting::Timers(Comms());
      buildInfo = new reporting::BuildInfo();
      state = new hemelb::lb::SimulationState(0.0001, 1000);
      net = new net::Net(Comms());
      latticeData = FourCubeLatticeData::Create(Comms(), 6, 5); // The 5 here is to match the topology size in the MPICommsMock
      LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(latticeData);
      latticeData->SwapOldAndNew(); //Needed since InitialiseAnisotropicTestData only initialises FOld
      cache = new lb::MacroscopicPropertyCache(*state, *latticeData);
      cache->densityCache.SetRefreshFlag();
      LbTestsHelper::UpdatePropertyCache<lb::lattices::D3Q15>(*latticeData, *cache, *state);
      incompChecker = new IncompressibilityCheckerMock(latticeData, net, state, *cache, *realTimers, 10.0);
      reporter = new Reporter("mock_path", "exampleinputfile");
      reporter->AddReportable(incompChecker);
      reporter->AddReportable(mockTimers);
      reporter->AddReportable(state);
      reporter->AddReportable(latticeData);
      reporter->AddReportable(buildInfo);

      auto AssertTemplate = [&](const std::string &expectation, const std::string &ttemplate) {
	ctemplate::StringToTemplateCache("TestFor" + ttemplate, ttemplate, ctemplate::DO_NOT_STRIP);
	std::string result;
	REQUIRE(ctemplate::ExpandTemplate("TestFor" + ttemplate,
					  ctemplate::DO_NOT_STRIP,
					  reporter->GetDictionary().GetRaw(),
					  &result));
	REQUIRE(expectation == result);
      };

      auto AssertValue = [&](const std::string & expectation, const std::string &symbol) {
	AssertTemplate(expectation, "{{" + symbol + "}}");
      };

      auto CheckTimingsTable = [&]() {
	std::stringstream expectation;
	expectation << std::setprecision(3);
	for (unsigned int row = 0; row < Timers::numberOfTimers; row++) {
	  expectation << "N" << TimersMock::timerNames[row] << "L" << row * 10.0 << "MI" << row * 15.0 << "ME"
		      << row * 2.0 << "MA" << row * 5.0 << " " << std::flush;
	}
	AssertTemplate(expectation.str(), "{{#TIMER}}N{{NAME}}L{{LOCAL}}MI{{MIN}}ME{{MEAN}}MA{{MAX}} {{/TIMER}}");
      };

      SECTION("TestInit") {
	AssertValue("exampleinputfile", "CONFIG");
      }

      

      SECTION("TestMainReport") {
	// Mock up some timings
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  for (unsigned int j = 0; j < i; j++) {
	    (*mockTimers)[i].Start();
	    (*mockTimers)[i].Stop();
	  }
	}
	mockTimers->Reduce(); // invoke the Timers MPI mock
	reporter->Image();
	reporter->Image();
	reporter->Image();
	for (unsigned int step = 0; step < 1000; step++) {
	  state->Increment();
	}
	REQUIRE(1001lu == state->GetTimeStep());
	reporter->FillDictionary();

	CheckTimingsTable();
	AssertTemplate("", "{{#UNSTABLE}} unstable{{/UNSTABLE}}");
	AssertTemplate("R0S64 R1S1000 R2S2000 R3S3000 R4S4000 ",
		       "{{#PROCESSOR}}R{{RANK}}S{{SITES}} {{/PROCESSOR}}");
	AssertTemplate(hemelb::reporting::mercurial_revision_number, "{{#BUILD}}{{REVISION}}{{/BUILD}}");
	AssertTemplate(hemelb::reporting::build_time, "{{#BUILD}}{{TIME}}{{/BUILD}}");
	AssertValue("3", "IMAGES");
	AssertValue("0.000100", "TIME_STEP_LENGTH");
	AssertValue("1000", "TOTAL_TIME_STEPS");
	AssertValue("1000", "STEPS");
	AssertValue("64", "SITES");
	AssertValue("1", "BLOCKS");
	AssertValue("216", "SITESPERBLOCK");
      }

      delete reporter;
      delete mockTimers;
      delete realTimers;
      delete cache;
      delete incompChecker;
      delete net;
      delete buildInfo;
    }
  }
}
