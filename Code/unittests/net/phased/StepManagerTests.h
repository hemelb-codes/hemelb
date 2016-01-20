
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "unittests/helpers/MockNetHelper.h"
#include "unittests/net/phased/MockConcern.h"
#include "unittests/net/phased/MockIteratedAction.h"
#include "net/phased/NetConcern.h"
#include <cppunit/TestFixture.h>

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      namespace phased
      {
        using namespace hemelb::net::phased;
        using namespace hemelb::unittests::helpers;
        class StepManagerTests : public CppUnit::TestFixture, public MockNetHelper
        {
            CPPUNIT_TEST_SUITE (StepManagerTests);
            CPPUNIT_TEST (TestConstruct);

            CPPUNIT_TEST (TestRegisterActor);

            CPPUNIT_TEST (TestRegisterAction);
            CPPUNIT_TEST (TestRegisterTwoConcerns);
            CPPUNIT_TEST (TestCallNonCommsActions);

            CPPUNIT_TEST (TestRegisterCommsConcern);
            CPPUNIT_TEST (TestCallCommsActions);

            CPPUNIT_TEST (TestRegisterThreeConcerns);
            CPPUNIT_TEST (TestCallPhaseActions);

            CPPUNIT_TEST (TestCallSpecialActions);
            CPPUNIT_TEST (TestCallAllActionsOnePhase);

            CPPUNIT_TEST (TestCallAllActionsManyPhases);
            CPPUNIT_TEST (TestCallAllActionsPhaseByPhase);

            CPPUNIT_TEST_SUITE_END();

          public:
            StepManagerTests() :
                MockNetHelper(),stepManager(NULL), action(NULL), concern(NULL), netConcern(NULL), action2(NULL), concern2(NULL)
            {
            }

            void setUp()
            {
              stepManager = new StepManager(3);

            }

            void tearDown()
            {
              MockNetHelper::tearDown();
              delete stepManager;
              delete netConcern;
              delete action;
              delete action2;
              delete concern;
              delete concern2;
            }

            void TestConstruct()
            {
              // PASS -- just assert the setup/teardown runs ok.
            }

            void TestRegisterActor()
            {
              action = new MockIteratedAction("mockOne");
              stepManager->RegisterIteratedActorSteps(*action);
              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 1u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 5u);
            }

            void TestRegisterAction()
            {
              concern = new MockConcern("mockTwo");
              stepManager->Register(0, steps::PreSend, *concern, 0);
              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 1u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 1u);
              stepManager->Register(0, steps::EndPhase, *concern, 1);
              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 1u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 2u);
            }

            void TestRegisterTwoConcerns()
            {
              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);

              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 2u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 7u);
            }

            void TestRegisterCommsConcern()
            {
              SetupMocks(4, 1);
              stepManager->RegisterCommsSteps(*netConcern, 0);
              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 1u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 3u);
            }

            void TestRegisterThreeConcerns()
            {
              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");
              SetupMocks(4, 1);

              stepManager->RegisterCommsSteps(*netConcern, 0);
              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);

              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 3u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 10u);
            }

            void TestCallCommsActions()
            {
              SetupMocks(4, 1);
              int payload0 = 0;
              int payload1 = 1;

              int payload0Expectation = 5;
              int payload1Expectation = 1;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              stepManager->RegisterCommsSteps(*netConcern, 0);

              stepManager->CallActionsForPhase(0);

              CPPUNIT_ASSERT_EQUAL(payload0, 5);
              netMock->ExpectationsAllCompleted();
            }

            void TestCallNonCommsActions()
            {
              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);

              stepManager->CallActionsForPhase(0);

              std::vector<int> shouldHaveCalled;

              shouldHaveCalled.push_back(0);
              shouldHaveCalled.push_back(1);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action->CallsSoFar());
            }

            void TestCallSpecialActions()
            {
              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::BeginAll, *concern, 37);

              stepManager->CallSpecialAction(steps::BeginAll);

              std::vector<int> shouldHaveCalled;
              shouldHaveCalled.push_back(37);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string(""), action->CallsSoFar());

              stepManager->CallSpecialAction(steps::EndAll);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("EndIteration, "), action->CallsSoFar());
            }

            void TestCallPhaseActions()
            {

              SetupMocks(4, 1);
              int payload0 = 0;
              int payload1 = 1;

              int payload0Expectation = 5;
              int payload1Expectation = 1;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);
              stepManager->RegisterCommsSteps(*netConcern, 0);

              stepManager->CallActionsForPhase(0);

              std::vector<int> shouldHaveCalled;

              shouldHaveCalled.push_back(0);
              shouldHaveCalled.push_back(1);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 5);
              netMock->ExpectationsAllCompleted();
            }

            void TestCallAllActionsOnePhase()
            {

              SetupMocks(4, 1);
              int payload0 = 0;
              int payload1 = 1;

              int payload0Expectation = 5;
              int payload1Expectation = 1;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              stepManager->RegisterIteratedActorSteps(*action);
              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);
              stepManager->Register(0, steps::BeginAll, *concern, 17);
              stepManager->RegisterCommsSteps(*netConcern, 0);

              stepManager->CallActions();

              std::vector<int> shouldHaveCalled;

              shouldHaveCalled.push_back(17);
              shouldHaveCalled.push_back(0);
              shouldHaveCalled.push_back(1);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 5);
              netMock->ExpectationsAllCompleted();
            }

            void TestCallAllActionsManyPhases()
            {

              SetupMocks(4, 1);
              int payload0 = 0;
              int payload1 = 1;

              int payload0Expectation = 5;
              int payload1Expectation = 1;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              action2 = new MockIteratedAction("mockThree");
              concern2 = new MockConcern("mockFour");

              stepManager->RegisterIteratedActorSteps(*action, 0);
              stepManager->RegisterIteratedActorSteps(*action2, 1);

              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);
              stepManager->Register(0, steps::BeginAll, *concern, 17);

              stepManager->Register(2, steps::PreSend, *concern2, 56);
              stepManager->Register(0, steps::EndPhase, *concern2, 42);
              stepManager->Register(0, steps::EndAll, *concern2, 13);

              stepManager->RegisterCommsForAllPhases(*netConcern);

              stepManager->CallActions();

              std::vector<int> shouldHaveCalled;

              shouldHaveCalled.push_back(17);
              shouldHaveCalled.push_back(0);
              shouldHaveCalled.push_back(1);

              std::vector<int> shouldHaveCalled2;

              shouldHaveCalled2.push_back(42); // gets called first cos is in phase 0
              shouldHaveCalled2.push_back(56);
              shouldHaveCalled2.push_back(13); // still called last cos is final special action

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   action2->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 5);
              netMock->ExpectationsAllCompleted();
            }

            void TestCallAllActionsPhaseByPhase()
            {

              SetupMocks(4, 1);
              int payload0 = 0;
              int payload1 = 1;

              action = new MockIteratedAction("mockOne");
              concern = new MockConcern("mockTwo");

              action2 = new MockIteratedAction("mockThree");
              concern2 = new MockConcern("mockFour");

              stepManager->RegisterIteratedActorSteps(*action, 0);
              stepManager->RegisterIteratedActorSteps(*action2, 1);

              stepManager->Register(0, steps::PreSend, *concern, 0);
              stepManager->Register(0, steps::EndPhase, *concern, 1);
              stepManager->Register(0, steps::BeginAll, *concern, 17);

              stepManager->Register(2, steps::PreSend, *concern2, 56);
              stepManager->Register(0, steps::EndPhase, *concern2, 42);
              stepManager->Register(0, steps::EndAll, *concern2, 13);

              stepManager->RegisterCommsForAllPhases(*netConcern);

              //-------------------BeginAll-----------------------------------------

              stepManager->CallSpecialAction(steps::BeginAll);

              std::vector<int> shouldHaveCalled;
              std::vector<int> shouldHaveCalled2;

              shouldHaveCalled.push_back(17);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string(""), action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string(""), action2->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 0);
              netMock->ExpectationsAllCompleted();

              //------------------Phase 0 ---------------------------------------------

              int payload0Expectation = 5;
              int payload1Expectation = 1;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              stepManager->CallActionsForPhase(0);

              shouldHaveCalled.push_back(0);
              shouldHaveCalled.push_back(1);

              shouldHaveCalled2.push_back(42); // gets called first cos is in phase 0

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string(""), action2->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 5);
              netMock->ExpectationsAllCompleted();

              // ------------------- Phase 1 ------------------------------------------

              payload0Expectation = 13;
              payload1Expectation = 4;
              payload1 = 4;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              stepManager->CallActionsForPhase(1);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action2->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 13);
              netMock->ExpectationsAllCompleted();

              // ------------------------- Phase 2 -----------------------------------------

              payload0Expectation = 77;
              payload1Expectation = 16;
              payload1 = 16;
              netMock->RequestSendR(payload1, 1);
              netMock->RequestReceiveR(payload0, 1);
              netMock->RequireSend(&payload1Expectation, 1, 1);
              netMock->RequireReceive(&payload0Expectation, 1, 1);

              stepManager->CallActionsForPhase(2);

              shouldHaveCalled2.push_back(56);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "),
                                   action2->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(payload0, 77);
              netMock->ExpectationsAllCompleted();

              // -------------------------- EndAll -----------------------------------------------

              stepManager->CallSpecialAction(steps::EndAll);

              shouldHaveCalled2.push_back(13); // still called last cos is final special action

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled2, concern2->ActionsCalled());

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   action->CallsSoFar());
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   action2->CallsSoFar());
              netMock->ExpectationsAllCompleted();
            }

            void SetupMocks(const proc_t core_count, const proc_t current_core)
            {
              MockNetHelper::setUp(core_count,current_core);
              netConcern = new NetConcern(*netMock);
            }

          private:

            StepManager *stepManager;
            MockIteratedAction *action;
            MockConcern *concern;
            NetConcern *netConcern;
            MockIteratedAction *action2;
            MockConcern *concern2;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION (StepManagerTests);
      }
    }
  }
}
#endif // ONCE
