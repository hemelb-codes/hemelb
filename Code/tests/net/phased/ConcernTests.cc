// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "net/phased/steps.h"

//#include "tests/helpers/Vectors.h"
#include "tests/net/phased/MockConcern.h"
#include "tests/net/phased/MockIteratedAction.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::net::phased;
    TEST_CASE("ConcernTests") {

      SECTION("TestCallAction") {
	MockConcern concern{"mockOne"};
	std::vector<int> shouldHaveCalled;

	concern.CallAction(0);
	shouldHaveCalled.push_back(0);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	concern.CallAction(1);
	shouldHaveCalled.push_back(1);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	concern.CallAction(3);
	shouldHaveCalled.push_back(3);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	concern.CallAction(1);
	shouldHaveCalled.push_back(1);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
      }

      SECTION("TestCallActorAsConcern") {
	MockIteratedAction actor{"mockTwo"};
	actor.CallAction(steps::BeginPhase);

	REQUIRE(std::string("RequestComms, ") == actor.CallsSoFar());

	actor.CallAction(steps::PreSend);
	REQUIRE(std::string("RequestComms, PreSend, ") == actor.CallsSoFar());

	actor.CallAction(steps::PreWait);
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, ") == actor.CallsSoFar());

	actor.CallAction(steps::EndPhase);
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== actor.CallsSoFar());

	actor.CallAction(steps::EndAll);
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== actor.CallsSoFar());
      }

    }
  }
}
