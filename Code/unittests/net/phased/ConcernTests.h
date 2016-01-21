
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_NET_PHASED_CONCERNTESTS_H
#define HEMELB_UNITTESTS_NET_PHASED_CONCERNTESTS_H

#include "unittests/helpers/CppUnitCompareVectors.h"
#include "unittests/net/phased/MockConcern.h"
#include "unittests/net/phased/MockIteratedAction.h"
#include "net/phased/steps.h"
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
            CPPUNIT_TEST (TestCallAction);
            CPPUNIT_TEST (TestCallActorAsConcern);
            CPPUNIT_TEST_SUITE_END();

          public:
            ConcernTests()
            {
            }

            void setUp()
            {
              concern = new MockConcern("mockOne");
              actor = new MockIteratedAction("mockTwo");
            }

            void tearDown()
            {
              delete concern;
              delete actor;
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

            void TestCallActorAsConcern()
            {
              actor->CallAction(steps::BeginPhase);

              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, "),actor->CallsSoFar());

              actor->CallAction(steps::PreSend);
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, "), actor->CallsSoFar());

              actor->CallAction(steps::PreWait);
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, "), actor->CallsSoFar());

              actor->CallAction(steps::EndPhase);
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, "), actor->CallsSoFar());

              actor->CallAction(steps::EndAll);
              CPPUNIT_ASSERT_EQUAL(std::string("RequestComms, PreSend, PreReceive, PostReceive, EndIteration, "),
                                   actor->CallsSoFar());
            }

          private:
            MockConcern * concern;
            MockIteratedAction * actor;

        };
        CPPUNIT_TEST_SUITE_REGISTRATION (ConcernTests);
      }
    }
  }
}
#endif // ONCE
