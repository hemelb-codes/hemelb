#ifndef HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "unittests/net/NetMock.h"
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
            CPPUNIT_TEST_SUITE_END();

          public:
            StepManagerTests() :
                netMock(NULL), communicatorMock(NULL), stepManager(NULL)
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
              concern = new MockIteratedAction("mockOne");
              stepManager->RegisterAllSteps(*concern);
              CPPUNIT_ASSERT_EQUAL(stepManager->ConcernCount(), 1u);
              CPPUNIT_ASSERT_EQUAL(stepManager->ActionCount(), 7u);
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
            MockIteratedAction *concern;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION (StepManagerTests);
      }
    }
  }
}
#endif // ONCE
