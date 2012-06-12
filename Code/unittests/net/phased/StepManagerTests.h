#ifndef HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "unittests/net/NetMock.h"
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
        class StepManagerTests : public CppUnit::TestFixture
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

            CPPUNIT_TEST_SUITE_END();

          public:
            StepManagerTests() :
                netMock(NULL), communicatorMock(NULL), stepManager(NULL), action(NULL), concern(NULL), netConcern(NULL)
            {
            }

            void setUp()
            {
              bool dummy;
              topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);
              stepManager = new StepManager();

            }

            void tearDown()
            {
              delete communicatorMock;
              delete netMock;
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
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 6u);
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
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 8u);
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
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 11u);
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
              stepManager->Register(0, steps::Reset, *concern, 19);

              stepManager->CallSpecialAction(steps::BeginAll);

              std::vector<int> shouldHaveCalled;
              shouldHaveCalled.push_back(37);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string(""), action->CallsSoFar());

              stepManager->CallSpecialAction(steps::Reset);
              shouldHaveCalled.push_back(19);

              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("Reset, "), action->CallsSoFar());

              stepManager->CallSpecialAction(steps::EndAll);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
              CPPUNIT_ASSERT_EQUAL(std::string("Reset, EndIteration, "), action->CallsSoFar());
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

            void SetupMocks(const proc_t core_count, const proc_t current_core)
            {
              communicatorMock = new topology::Communicator(current_core, core_count);
              netMock = new net::NetMock(*communicatorMock);
              netConcern = new NetConcern(*netMock);
            }

          private:

            net::NetMock *netMock;
            topology::Communicator *communicatorMock;
            StepManager *stepManager;
            MockIteratedAction *action;
            MockConcern *concern;
            NetConcern *netConcern;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION (StepManagerTests);
      }
    }
  }
}
#endif // ONCE
