// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "configuration/CommandLine.h"
#include "io/PathManager.h"

#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace tests
  {

    std::unique_ptr<io::PathManager> setup(const char* confXml) {
      const int processorCount = 5;
      const int argc = 7;
      const char* argv[] = {
	"hemelb",
	"-in",
	confXml,
	"-i",
	"1",
	"-ss",
	"1111",
      };

      auto cl = configuration::CommandLine(argc, argv);
      return std::make_unique<io::PathManager>(cl, true, processorCount);
    }

    TEST_CASE_METHOD(helpers::FolderTestFixture, "PathManager") {
      SECTION("Local config XML path") {
	auto fileManager = setup("config.xml");

	SECTION("Create") {
	  AssertPresent("results");
	  AssertPresent("results/Images");
	}

	SECTION("NameInvention") {
	  REQUIRE(std::string("./results") == fileManager->GetReportPath());
	}

      }

      SECTION("Abs path to config XML") {
	const std::string targetConfig = GetTempdir() + "/config.xml"; // note this resource doesn't exist -- not a problem
	ReturnToOrigin(); // even if we're not in current dir, explicit path should cause the results to be created in the tmpdir
	auto fileManager = setup(targetConfig.c_str());
	MoveToTempdir(); // go back to the tempdir and check the files were created in the right place

	SECTION("Create") {
	  AssertPresent("results");
	  AssertPresent("results/Images");
	}

	SECTION("NameInvention") {
	  REQUIRE(GetTempdir() + "/results" == fileManager->GetReportPath());
	}
      }

    }

  }
}
