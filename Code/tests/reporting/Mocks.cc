// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/reporting/Mocks.h"

#include <catch2/catch.hpp>

#include "net/IOCommunicator.h"
#include "reporting/Timers.h"

namespace hemelb::tests
{
    double ClockMock::operator()() {
        fakeTime += 10.0;
        return fakeTime;
    }

    std::vector<double> MPICommsMock::Reduce(const std::vector<double> & sendbuf, MPI_Op, int root) {
        auto const count = sendbuf.size();
        REQUIRE(reporting::Timers::last == count);
        std::vector<double> recvbuf(count);

        for (int i = 0; i < count; i++) {
            REQUIRE(10.0 * i == sendbuf[i]);
            recvbuf[i] = 5.0 * i * calls;
        }
        calls++;
        return recvbuf;
    }

    int MPICommsMock::Size() const
    {
        return 5;
    }

}

