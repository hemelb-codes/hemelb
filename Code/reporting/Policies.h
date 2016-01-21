
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_POLICIES_H
#define HEMELB_REPORTING_POLICIES_H

/**
 * @file
 * This file contains Policy classes defining how the report generator should interact with the file system, MPI, and the clock.
 * These are used as template arguments to the principal classes.
 * In file unittests/reporting/Mocks.h, mock versions of these policies are defined to facilitate testing.
 */

#include <fstream>
#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace reporting
  {

    /**
     * Abstraction of interaction with MPI system.
     * Mocked by hemelb::unittests::reporting::MPICommsMock
     */
    class MPICommsPolicy
    {
      public:
        /**
         * Stores a pointer to the MPI Network topology singleton.
         */
        MPICommsPolicy(const net::IOCommunicator& inst_) :
            instance(inst_)
        {
        }
      protected:
        /**
         * Wrap the MPI reduce call.
         * @param sendbuf Pointer to start of send buffer
         * @param recvbuf Pointer to start of receive buffer
         * @param count Number of data to send/receive
         * @param datatype MPI Datatype
         * @param op Operator to reduce under
         * @param root Target processor where the reduction is reduced to.
         * @param comm MPI Communicator to reduce over
         * @return
         */
        int Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root)
        {
          return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, instance);
        }

        /**
         * Total number of MPI nodes in the communicator.
         * @return Total number of processors, some of which may be shared on a single machine.
         */
        proc_t GetProcessorCount() const
        {
          return instance.Size();
        }

      private:
        const net::IOCommunicator& instance; //! Reference to the singleton instance of the MPI topology
    };

    /**
     * A way to get the time.
     * Mocked by hemelb::unittests::reporting::ClockMock
     */
    class HemeLBClockPolicy
    {
      protected:
        /**
         * Get the time
         * @return current time in seconds.
         */
        static double CurrentTime()
        {
          return hemelb::util::myClock();
        }
    };
  }
}
#endif // ONCE
