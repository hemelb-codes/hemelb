#ifndef HEMELB_REPORTING_POLICIES_H
#define HEMELB_REPORTING_POLICIES_H

/***
 * This file contains Policy classes defining how the report generator should interact with the file system, MPI, and the clock.
 * These are used as template arguments to the principal classes.
 * In file unittests/reporting/Mocks.h, mock versions of these policies are defined to facilitate testing.
 */

#include <stdarg.h>
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace reporting
  {
    /***
     * Policy defining writing text output to a simple C-style file.
     * Mocked by hemelb::unittests::reporting::WriterMock
     */
    class FileWriterPolicy
    {
      public:
        /**
         * Opens a file at the specified path.
         * @param path Full or relative path to where the file should be opened.
         */
        FileWriterPolicy(const std::string &path):
          file(fopen(path.c_str(), "w"))
        {
        }
        /**
         * Closes the open file.
         */
        ~FileWriterPolicy()
        {
          fclose(file);
        }
      protected:
        /**
         * Wraps printf-style behaviour.
         * Uses variable argument list and a format string to write to the file.
         * @param format C printf style format string.
         */
        void Print(const char * format, ...) const
        {
          std::va_list arg;
          va_start(arg, format);
          vfprintf(file, format, arg);
          va_end(arg);
        }
      private:
        FILE* const file; //! Handle of file used.
    };

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
        MPICommsPolicy() :
            instance(*hemelb::topology::NetworkTopology::Instance())
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
        int Reduce(void *sendbuf,
                   void *recvbuf,
                   int count,
                   MPI_Datatype datatype,
                   MPI_Op op,
                   int root,
                   MPI_Comm comm)
        {
          return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        }
        /**
         * Number of fluid sites assigned to processor.
         * @param n Processor number to be investigated.
         * @return Number of fluid sites assigned to that processor.
         */
        size_t FluidSitesOnProcessor(int n) const
        {
          return instance.FluidSitesOnEachProcessor[n];
        }
        /**
         * Total number of MPI nodes in the communicator.
         * @return Total number of processors, some of which may be shared on a single machine.
         */
        proc_t GetProcessorCount() const
        {
          return instance.GetProcessorCount();
        }
        /**
         * Number of machines, as opposed to processors.
         * @return Number of machines.
         */
        unsigned int GetMachineCount() const
        {
          return instance.GetMachineCount();
        }
        /**
         * Max depth of all processors within the machine topology.
         * @return Greatest processor depth.
         */
        int GetDepths() const
        {
          return instance.GetDepths();
        }

      private:
        const hemelb::topology::NetworkTopology& instance; //! Reference to the singleton instance of the MPI topology
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
