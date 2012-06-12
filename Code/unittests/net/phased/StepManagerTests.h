#ifndef HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "unittests/net/NetMock.h"
#include "unittests/net/phased/MockConcern.h"
#include "unittests/net/phased/MockIteratedAction.h"
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
            CPPUNIT_TEST (TestCallActions);
            CPPUNIT_TEST_SUITE_END();

          public:
            StepManagerTests() :
                netMock(NULL), communicatorMock(NULL), stepManager(NULL), action(NULL), concern(NULL)
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

            void TestCallActions()
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

            void SetupMocks(const proc_t core_count, const proc_t current_core)
            {
              communicatorMock = new topology::Communicator(current_core, core_count);
              netMock = new net::NetMock(*communicatorMock);

            }

          private:

            net::NetMock *netMock;
            topology::Communicator *communicatorMock;
            StepManager *stepManager;
            MockIteratedAction *action;
            MockConcern *concern;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION (StepManagerTests);
      }
    }
  }
}
#endif // ONCE
