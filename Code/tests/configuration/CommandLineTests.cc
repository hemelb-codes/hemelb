// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>
#include <catch2/catch.hpp>

#include "configuration/CommandLine.h"
#include "resources/Resource.h"

#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace helpers;

    TEST_CASE_METHOD(helpers::FolderTestFixture, "CommandLine") {
      auto configFile = resources::Resource("four_cube.xml").Path();
      const int argc = 7;
      const char* argv[] = {
	"hemelb",
	"-in",
	configFile.c_str(),
	"-i",
	"1",
	"-ss",
	"1111"
      };
      
        
      SECTION("Construct"){
	auto options = std::make_unique<hemelb::configuration::CommandLine>(argc, argv);
	REQUIRE(options != nullptr);
      }
    }
  }
}
