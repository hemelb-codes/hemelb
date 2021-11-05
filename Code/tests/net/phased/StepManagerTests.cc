// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "net/phased/StepManager.h"
#include "net/phased/NetConcern.h"

#include "tests/helpers/MockNetHelper.h"
#include "tests/net/phased/MockConcern.h"
#include "tests/net/phased/MockIteratedAction.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::net::phased;
    using namespace hemelb::tests::helpers;

    TEST_CASE_METHOD(MockNetHelper, "StepManagerTests") {
      StepManager stepManager{3};
      // MockIteratedAction *action = nullptr;
      // MockConcern *concern = nullptr;
      // NetConcern *netConcern = nullptr;
      // MockIteratedAction *action2 = nullptr;
      // MockConcern *concern2 = nullptr;
      const proc_t core_count = 4;
      const proc_t current_core = 1;
      setUp(core_count, current_core);
      NetConcern netConcern{*netMock};
      
      SECTION("TestRegisterActor") {
	MockIteratedAction action{"mockOne"};
      	stepManager.RegisterIteratedActorSteps(action);
      	REQUIRE(stepManager.ConcernCount() == 1u);
      	REQUIRE(stepManager.ActionCount() == 5u);
      }

      SECTION("TestRegisterAction") {
      	MockConcern concern("mockTwo");
      	stepManager.Register(0, steps::PreSend, concern, 0);
	REQUIRE(stepManager.ConcernCount() == 1u);
	REQUIRE(stepManager.ActionCount() == 1u);
      	stepManager.Register(0, steps::EndPhase, concern, 1);
	REQUIRE(stepManager.ConcernCount() == 1u);
	REQUIRE(stepManager.ActionCount() == 2u);
      }

      SECTION("TestRegisterTwoConcerns") {
      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);

	REQUIRE(stepManager.ConcernCount() == 2u);
	REQUIRE(stepManager.ActionCount() == 7u);
      }

      SECTION("TestRegisterCommsConcern") {
      	stepManager.RegisterCommsSteps(netConcern, 0);
	REQUIRE(stepManager.ConcernCount() == 1u);
	REQUIRE(stepManager.ActionCount() == 3u);
      }

      SECTION("TestRegisterThreeConcerns") {
      	MockIteratedAction action("mockOne");
	MockConcern concern("mockTwo");
      	
      	stepManager.RegisterCommsSteps(netConcern, 0);
      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);

	REQUIRE(stepManager.ConcernCount() == 3u);
	REQUIRE(stepManager.ActionCount() == 10u);
      }

      SECTION("TestCallCommsActions") {
      	int payload0 = 0;
      	int payload1 = 1;

      	int payload0Expectation = 5;
      	int payload1Expectation = 1;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	stepManager.RegisterCommsSteps(netConcern, 0);

      	stepManager.CallActionsForPhase(0);

	REQUIRE(payload0 == 5);
      	netMock->ExpectationsAllCompleted();
      }

      SECTION("TestCallNonCommsActions") {
	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);

      	stepManager.CallActionsForPhase(0);

      	std::vector<int> shouldHaveCalled;

      	shouldHaveCalled.push_back(0);
      	shouldHaveCalled.push_back(1);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action.CallsSoFar());
      }

      SECTION("TestCallSpecialActions") {
      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::BeginAll, concern, 37);

      	stepManager.CallSpecialAction(steps::BeginAll);

      	std::vector<int> shouldHaveCalled;
      	shouldHaveCalled.push_back(37);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(std::string("") == action.CallsSoFar());

      	stepManager.CallSpecialAction(steps::EndAll);
	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(std::string("EndIteration, ") == action.CallsSoFar());
      }

      SECTION("TestCallPhaseActions") {
      	int payload0 = 0;
      	int payload1 = 1;

      	int payload0Expectation = 5;
      	int payload1Expectation = 1;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);
      	stepManager.RegisterCommsSteps(netConcern, 0);

      	stepManager.CallActionsForPhase(0);

      	std::vector<int> shouldHaveCalled;

      	shouldHaveCalled.push_back(0);
      	shouldHaveCalled.push_back(1);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action.CallsSoFar());
	REQUIRE(payload0 == 5);
      	netMock->ExpectationsAllCompleted();
      }

      SECTION("TestCallAllActionsOnePhase") {
      	int payload0 = 0;
      	int payload1 = 1;

      	int payload0Expectation = 5;
      	int payload1Expectation = 1;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	stepManager.RegisterIteratedActorSteps(action);
      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);
      	stepManager.Register(0, steps::BeginAll, concern, 17);
      	stepManager.RegisterCommsSteps(netConcern, 0);

      	stepManager.CallActions();

      	std::vector<int> shouldHaveCalled;

      	shouldHaveCalled.push_back(17);
      	shouldHaveCalled.push_back(0);
      	shouldHaveCalled.push_back(1);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());

	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== action.CallsSoFar());
      	REQUIRE(payload0 == 5);
      	netMock->ExpectationsAllCompleted();
      }

      SECTION("TestCallAllActionsManyPhases") {
      	int payload0 = 0;
      	int payload1 = 1;

      	int payload0Expectation = 5;
      	int payload1Expectation = 1;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	MockIteratedAction action2("mockThree");
      	MockConcern concern2("mockFour");

      	stepManager.RegisterIteratedActorSteps(action, 0);
      	stepManager.RegisterIteratedActorSteps(action2, 1);

      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);
      	stepManager.Register(0, steps::BeginAll, concern, 17);

      	stepManager.Register(2, steps::PreSend, concern2, 56);
      	stepManager.Register(0, steps::EndPhase, concern2, 42);
      	stepManager.Register(0, steps::EndAll, concern2, 13);

      	stepManager.RegisterCommsForAllPhases(netConcern);

      	stepManager.CallActions();

      	std::vector<int> shouldHaveCalled;

      	shouldHaveCalled.push_back(17);
      	shouldHaveCalled.push_back(0);
      	shouldHaveCalled.push_back(1);

      	std::vector<int> shouldHaveCalled2;

      	shouldHaveCalled2.push_back(42); // gets called first cos is in phase 0
      	shouldHaveCalled2.push_back(56);
      	shouldHaveCalled2.push_back(13); // still called last cos is final special action

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());

	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== action.CallsSoFar());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== action2.CallsSoFar());
	REQUIRE(payload0 == 5);
      	netMock->ExpectationsAllCompleted();
      }

      SECTION("TestCallAllActionsPhaseByPhase") {
	int payload0 = 0;
      	int payload1 = 1;

      	MockIteratedAction action("mockOne");
      	MockConcern concern("mockTwo");

      	MockIteratedAction action2("mockThree");
      	MockConcern concern2("mockFour");

      	stepManager.RegisterIteratedActorSteps(action, 0);
      	stepManager.RegisterIteratedActorSteps(action2, 1);

      	stepManager.Register(0, steps::PreSend, concern, 0);
      	stepManager.Register(0, steps::EndPhase, concern, 1);
      	stepManager.Register(0, steps::BeginAll, concern, 17);

      	stepManager.Register(2, steps::PreSend, concern2, 56);
      	stepManager.Register(0, steps::EndPhase, concern2, 42);
      	stepManager.Register(0, steps::EndAll, concern2, 13);

      	stepManager.RegisterCommsForAllPhases(netConcern);

      	//-------------------BeginAll-----------------------------------------

      	stepManager.CallSpecialAction(steps::BeginAll);

      	std::vector<int> shouldHaveCalled;
      	std::vector<int> shouldHaveCalled2;

      	shouldHaveCalled.push_back(17);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());
	REQUIRE(std::string("") == action.CallsSoFar());
	REQUIRE(std::string("") == action2.CallsSoFar());
	REQUIRE(payload0 == 0);
      	netMock->ExpectationsAllCompleted();

      	//------------------Phase 0 ---------------------------------------------

      	int payload0Expectation = 5;
      	int payload1Expectation = 1;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	stepManager.CallActionsForPhase(0);

      	shouldHaveCalled.push_back(0);
      	shouldHaveCalled.push_back(1);

      	shouldHaveCalled2.push_back(42); // gets called first cos is in phase 0

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action.CallsSoFar());
	REQUIRE(std::string("") == action2.CallsSoFar());
	REQUIRE(payload0 == 5);
      	netMock->ExpectationsAllCompleted();

      	// ------------------- Phase 1 ------------------------------------------

      	payload0Expectation = 13;
      	payload1Expectation = 4;
      	payload1 = 4;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	stepManager.CallActionsForPhase(1);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action.CallsSoFar());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action2.CallsSoFar());
	REQUIRE(payload0 == 13);
      	netMock->ExpectationsAllCompleted();

      	// ------------------------- Phase 2 -----------------------------------------

      	payload0Expectation = 77;
      	payload1Expectation = 16;
      	payload1 = 16;
      	netMock->RequestSendR(payload1, 1);
      	netMock->RequestReceiveR(payload0, 1);
      	netMock->RequireSend(&payload1Expectation, 1, 1);
      	netMock->RequireReceive(&payload0Expectation, 1, 1);

      	stepManager.CallActionsForPhase(2);

      	shouldHaveCalled2.push_back(56);

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action.CallsSoFar());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, ")
		== action2.CallsSoFar());
	REQUIRE(payload0 == 77);
      	netMock->ExpectationsAllCompleted();

      	// -------------------------- EndAll -----------------------------------------------

      	stepManager.CallSpecialAction(steps::EndAll);

      	shouldHaveCalled2.push_back(13); // still called last cos is final special action

	REQUIRE(shouldHaveCalled == concern.ActionsCalled());
	REQUIRE(shouldHaveCalled2 == concern2.ActionsCalled());

	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== action.CallsSoFar());
	REQUIRE(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, ")
		== action2.CallsSoFar());
      	netMock->ExpectationsAllCompleted();
      }
    }

  }
}
