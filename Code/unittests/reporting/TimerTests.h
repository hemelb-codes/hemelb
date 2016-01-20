
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_REPORTING_TIMERTESTS_H
#define HEMELB_UNITTESTS_REPORTING_TIMERTESTS_H

#include <cppunit/TestFixture.h>
#include "reporting/Timers.h"
#include "reporting/Timers.hpp"
#include "unittests/reporting/Mocks.h"
#include "unittests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace reporting
    {
      using namespace hemelb::reporting;
      class TimerTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE(TimerTests);
          CPPUNIT_TEST(TestInitialization);
          CPPUNIT_TEST(TestStartStop);
          CPPUNIT_TEST(TestSetTime);
          CPPUNIT_TEST(TestMultipleStartStop);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            timer = new TimerBase<ClockMock>();
          }

          void tearDown()
          {
            delete timer;
          }

          void TestInitialization()
          {
            // should default to zero time on startup
            CPPUNIT_ASSERT_EQUAL(0.0, timer->Get());
          }
          void TestStartStop()
          {
            timer->Start(); // clock mock at 10.0
            timer->Stop(); // clock mock at 20.0
            CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0, timer->Get(), 1e-6);
          }
          void TestMultipleStartStop()
          {
            timer->Start(); // clock mock at 10.0
            timer->Stop(); // clock mock at 20.0
            timer->Start(); // clock mock at 30.0
            timer->Stop(); // clock mock at 40.0
            CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0, timer->Get(), 1e-6);
          }
          void TestSetTime()
          {
            timer->Set(15.0);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0, timer->Get(), 1e-6);
            timer->Start(); // clock mock at 10.0
            timer->Stop(); // clock mock at 20.0
            CPPUNIT_ASSERT_DOUBLES_EQUAL(25.0, timer->Get(), 1e-6);
          }

        private:
          TimerBase<ClockMock> *timer;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(TimerTests);

      class TimersTests : public helpers::HasCommsTestFixture
      {
          CPPUNIT_TEST_SUITE(TimersTests);
          CPPUNIT_TEST(TestInitialization);
          CPPUNIT_TEST(TestTimersSeparate);
          CPPUNIT_TEST(TestReduce);
          CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {
            timers = new TimersBase<ClockMock, MPICommsMock>(Comms());
          }

          void tearDown()
          {
            delete timers;
            helpers::HasCommsTestFixture::tearDown();
          }

          void TestInitialization()
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (*timers)[Timers::total].Get(), 1e-6);
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, (*timers)[i].Get(), 1e-6);
            }
          }

          void TestTimersSeparate()
          {
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              for (unsigned int j = 0; j < i; j++)
              {
                (*timers)[i].Start();
                (*timers)[i].Stop();
              }
            }
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0 * i, (*timers)[i].Get(), 1e-6);
            }
          }

          void TestReduce()
          {
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              for (unsigned int j = 0; j < i; j++)
              {
                (*timers)[i].Start();
                (*timers)[i].Stop();
              }
            }
            timers->Reduce(); // trigger the mock
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(i * 5.0, timers->Maxes()[i], 1e-6);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(i * 2.0, timers->Means()[i], 1e-6);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(i * 15.0, timers->Mins()[i], 1e-6);
            }
          }

        private:
          TimersBase<ClockMock, MPICommsMock> *timers;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION(TimersTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_REPORTING_TIMERTESTS_H_ */
