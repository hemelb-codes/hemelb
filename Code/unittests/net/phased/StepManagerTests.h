#ifndef HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_STEPMANAGERTESTS_H
#include "net/phased/StepManager.h"
#include "unittests/net/NetMock.h"
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
            CPPUNIT_TEST (TestConstruct);CPPUNIT_TEST_SUITE_END();

          public:
            StepManagerTests() :
                netMock(NULL), communicatorMock(NULL)
            {
            }

            void setUp()
            {
              bool dummy;
              topology::NetworkTopology::Instance()->Init(0, NULL, &dummy);

            }

            void tearDown()
            {
              delete communicatorMock;
              delete netMock;
            }

            void TestConstruct()
            {
            }

            void SetupMocks(const proc_t core_count, const proc_t current_core)
            {
              communicatorMock = new topology::Communicator(current_core, core_count);
              netMock = new net::NetMock(*communicatorMock);

            }

          private:

            net::NetMock *netMock;
            topology::Communicator *communicatorMock;
        };

        CPPUNIT_TEST_SUITE_REGISTRATION (StepManagerTests);
      }
    }
  }
}
#endif // ONCE
