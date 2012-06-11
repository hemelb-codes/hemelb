#ifndef HEMELB_UNITTESTS_NET_PHASED_CONCERNTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_CONCERNSTESTS_H
#include "net/phased/Concern.h"
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
        class ConcernTests : public CppUnit::TestFixture
        {
            CPPUNIT_TEST_SUITE (ConcernTests);
            CPPUNIT_TEST (TestConstruct);CPPUNIT_TEST_SUITE_END();

          public:
            ConcernTests()
            {
            }

            void setUp()
            {

            }

            void tearDown()
            {
            }

            void TestConstruct()
            {
            }


          private:

        };

        CPPUNIT_TEST_SUITE_REGISTRATION (ConcernTests);
      }
    }
  }
}
#endif // ONCE
