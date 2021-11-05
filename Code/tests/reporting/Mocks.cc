// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/reporting/Mocks.h"

#include <catch2/catch.hpp>

#include "net/IOCommunicator.h"
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

    MPICommsMock::MPICommsMock(const net::IOCommunicator& ignored) :
              calls(1)
    {
    }
    int MPICommsMock::Reduce(double *sendbuf,
			     double *recvbuf,
			     int count,
			     MPI_Datatype datatype,
			     MPI_Op op,
			     int root) {
      REQUIRE((int)reporting::Timers::last == count);
      for (int i = 0; i < count; i++) {
	REQUIRE(10.0 * i == sendbuf[i]);
	recvbuf[i] = 5.0 * i * calls;
      }
      calls++;
      return 0;
    }
    proc_t MPICommsMock::GetProcessorCount()
    {
      return 5;
    }

  }
}

