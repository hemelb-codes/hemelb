#ifndef HEMELB_UNITTESTS_NET_PHASED_CONCERNTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_CONCERNTESTS_H

#include "unittests/helpers/CppUnitCompareVectors.h"
#include "unittests/net/phased/MockConcern.h"
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
            CPPUNIT_TEST (TestConstruct);
            CPPUNIT_TEST (TestCallAction);CPPUNIT_TEST_SUITE_END();

          public:
            ConcernTests()
            {
            }

            void setUp()
            {
              concern = new MockConcern();
            }

            void tearDown()
            {
            }

            void TestConstruct()
            {
              // PASS -- just verify setUp and tearDown
            }

            void TestCallAction()
            {
              std::vector<int> shouldHaveCalled;

              concern->CallAction(0);
              shouldHaveCalled.push_back(0);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              concern->CallAction(1);
              shouldHaveCalled.push_back(1);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              concern->CallAction(3);
              shouldHaveCalled.push_back(3);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());

              concern->CallAction(1);
              shouldHaveCalled.push_back(1);
              CPPUNIT_ASSERT_EQUAL(shouldHaveCalled, concern->ActionsCalled());
            }

          private:
            MockConcern * concern;

        };
        CPPUNIT_TEST_SUITE_REGISTRATION (ConcernTests);
      }
    }
  }
}
#endif // ONCE
