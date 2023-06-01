// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REPORTING_MOCKS_H
#define HEMELB_TESTS_REPORTING_MOCKS_H

#include <mpi.h>
#include "units.h"

namespace hemelb::net {
    class IOCommunicator;
}

namespace hemelb::tests
{
    class ClockMock {
    public:
      ClockMock() = default;
      double operator()();
    private:
      double fakeTime = 0.0;
    };

    class MPICommsMock {
        unsigned calls = 1;
    public:
        int Size() const;
        std::vector<double> Reduce(std::vector<double> const&, MPI_Op, int root);
    };
}

#endif // HEMELB_TESTS_REPORTING_MOCKS_H
