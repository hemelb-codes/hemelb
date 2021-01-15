// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "resources/Resource.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/LaddFail.h"

namespace hemelb
{
  namespace tests
  {
    using namespace configuration;
    namespace
    {
      // For section XMLFileContent
      struct CfgChecker {
	using result_type = bool;
	template<typename T>
	bool operator()(T) const {
	  return false;
	}
	bool operator()(const EquilibriumIC& eqIC) const {
	  return 80.0 == eqIC.p_mmHg;
	}
      };
    }
    TEST_CASE_METHOD(helpers::FolderTestFixture, "SimConfig") {

      SECTION("0_2_0_Read") {
	LADD_FAIL();
	// smoke test the configuration as having loaded OK
	auto config = std::unique_ptr<SimConfig>(SimConfig::New(resources::Resource("config0_2_0.xml").Path()));
	REQUIRE(3000lu == config->GetTotalTimeSteps());
	REQUIRE(0.0001 == config->GetTimeStepLength());
	lb::iolets::InOutLetCosine* inlet = dynamic_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0]);
	REQUIRE(inlet != nullptr);
	REQUIRE(Approx(6000.0) == inlet->GetPeriod());

	// Check that in the absence of the <monitoring> XML element things get initiliased properly
	const hemelb::configuration::MonitoringConfig* monConfig = config->GetMonitoringConfiguration();
	REQUIRE(!monConfig->doConvergenceCheck);
	REQUIRE(!monConfig->doIncompressibilityCheck);
	REQUIRE(!monConfig->convergenceTerminate);
	REQUIRE(Approx(0.0) == monConfig->convergenceRelativeTolerance);
      }

      SECTION("0_2_1_Read") {
	LADD_FAIL();
	// smoke test the configuration as having loaded OK
	auto config = std::unique_ptr<SimConfig>(SimConfig::New(resources::Resource("config.xml").Path()));
	REQUIRE(3000lu == config->GetTotalTimeSteps());
	REQUIRE(0.0001 == config->GetTimeStepLength());
	lb::iolets::InOutLetCosine* inlet = dynamic_cast<lb::iolets::InOutLetCosine*>(config->GetInlets()[0]);
	REQUIRE(inlet != nullptr);
	REQUIRE(Approx(6000.0) == inlet->GetPeriod());

	const hemelb::configuration::MonitoringConfig* monConfig = config->GetMonitoringConfiguration();
	REQUIRE(monConfig->doConvergenceCheck);
	REQUIRE(monConfig->doIncompressibilityCheck);
	REQUIRE(monConfig->convergenceTerminate);
	REQUIRE(1e-9 == monConfig->convergenceRelativeTolerance);
	REQUIRE(monConfig->convergenceVariable == extraction::OutputField::Velocity);
	REQUIRE(0.01 == monConfig->convergenceReferenceValue); // 1 m/s * (delta_t / delta_x) = 0.01
      }
      
      SECTION("XMLFileContent") {
	LADD_FAIL();
	//Round trip the config twice.
	CopyResourceToTempdir("config.xml");
	auto config = std::unique_ptr<SimConfig>(SimConfig::New("config.xml"));

	auto& ICconfig = config->GetInitialCondition();
	REQUIRE(boost::apply_visitor(CfgChecker{}, ICconfig));
      }

    }
  }
}
