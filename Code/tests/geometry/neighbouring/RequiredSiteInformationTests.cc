// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "geometry/neighbouring/RequiredSiteInformation.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::geometry::neighbouring;
    TEST_CASE("RequiredSiteInformationTests") {
      RequiredSiteInformation requirements;

      SECTION("TestConstruct")
	{
	  REQUIRE(!requirements.IsRequired(terms::Distribution));
	}

      SECTION("TestSpecifyRequirement")
	{
	  requirements.Require(terms::Distribution);
	  REQUIRE(requirements.IsRequired(terms::Distribution));
	}

      SECTION("TestAny")
	{
	  REQUIRE(!requirements.RequiresAny());
	  requirements.Require(terms::Distribution);
	  REQUIRE(requirements.RequiresAny());
	}

      SECTION("TestAnyNonFieldDependent")
	{
	  REQUIRE(!requirements.RequiresAnyNonFieldDependent());
	  requirements.Require(terms::Distribution);
	  REQUIRE(!requirements.RequiresAnyNonFieldDependent());
	  requirements.Require(terms::SiteData);
	  REQUIRE(requirements.RequiresAnyNonFieldDependent());
	}

      SECTION("TestAnyFieldDependent")
	{
	  REQUIRE(!requirements.RequiresAnyFieldDependent());
	  requirements.Require(terms::SiteData);
	  REQUIRE(!requirements.RequiresAnyFieldDependent());
	  requirements.Require(terms::Distribution);
	  REQUIRE(requirements.RequiresAnyFieldDependent());
	}

      SECTION("TestAnyMacroscopic")
	{
	  REQUIRE(!requirements.RequiresAnyMacroscopic());
	  requirements.Require(terms::SiteData);
	  REQUIRE(!requirements.RequiresAnyMacroscopic());
	  requirements.Require(terms::Distribution);
	  REQUIRE(!requirements.RequiresAnyMacroscopic());
	  requirements.Require(terms::Velocity);
	  REQUIRE(requirements.RequiresAnyMacroscopic());
	}
    }
  }
}
