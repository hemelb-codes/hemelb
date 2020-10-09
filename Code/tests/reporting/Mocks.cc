
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/reporting/Mocks.h"

#include <catch2/catch.hpp>

#include "reporting/Timers.h"

namespace hemelb
{
  namespace tests
  {
    ClockMock::ClockMock() : fakeTime(0)
    {
    }
    double ClockMock::CurrentTime() {
      fakeTime += 10.0;
      return fakeTime;
    }

  }
}

