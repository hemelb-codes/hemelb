// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "reporting/Timers.h"
#include "reporting/Timers.hpp"

#include "tests/reporting/Mocks.h"
#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::reporting;

    TEST_CASE("TimerTests") {
      auto timer = TimerBase<ClockMock>();

      SECTION("TestInitialization") {
	// should default to zero time on startup
	REQUIRE(0.0 == timer.Get());
      }

      SECTION("TestStartStop") {
	timer.Start(); // clock mock at 10.0
	timer.Stop(); // clock mock at 20.0
	REQUIRE(Approx(10.0) == timer.Get());
      }

      SECTION("TestMultipleStartStop") {
	timer.Start(); // clock mock at 10.0
	timer.Stop(); // clock mock at 20.0
	timer.Start(); // clock mock at 30.0
	timer.Stop(); // clock mock at 40.0
	REQUIRE(Approx(20.0) == timer.Get());
      }

      SECTION("TestSetTime") {
	timer.Set(15.0);
	REQUIRE(Approx(15.0) == timer.Get());
	timer.Start(); // clock mock at 10.0
	timer.Stop(); // clock mock at 20.0
	REQUIRE(Approx(25.0) == timer.Get());
      }
    }

    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "TimersTests") {
      auto timers = TimersBase<ClockMock, MPICommsMock>(Comms());
      
      SECTION("TestInitialization") {
	REQUIRE(Approx(0.0) == timers[Timers::total].Get());
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  REQUIRE(Approx(0.0) == timers[i].Get());
	}
      }

      SECTION("TestTimersSeparate") {
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  for (unsigned int j = 0; j < i; j++) {
	    timers[i].Start();
	    timers[i].Stop();
	  }
	}
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  REQUIRE(Approx(10.0 * i) == timers[i].Get());
	}
      }

      SECTION("TestReduce") {
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  for (unsigned int j = 0; j < i; j++) {
	    timers[i].Start();
	    timers[i].Stop();
	  }
	}
	timers.Reduce(); // trigger the mock
	for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
	  REQUIRE(Approx(i * 5.0) == timers.Maxes()[i]);
	  REQUIRE(Approx(i * 2.0) == timers.Means()[i]);
	  REQUIRE(Approx(i * 15.0) == timers.Mins()[i]);
	}
      }

    }

  }
}
