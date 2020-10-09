
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REPORTING_MOCKS_H
#define HEMELB_TESTS_REPORTING_MOCKS_H

#include <mpi.h>
#include "units.h"

namespace hemelb
{
  namespace tests
  {
    class ClockMock {
    public:
      ClockMock();
    protected:
      double CurrentTime();
    private:
      double fakeTime;
    };
  }
}

#endif // HEMELB_TESTS_REPORTING_MOCKS_H
