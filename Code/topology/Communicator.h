#ifndef HEMELB_TOPOLOGY_COMMUNICATOR_H
#define HEMELB_TOPOLOGY_COMMUNICATOR_H

#include "units.h"
#include <mpi.h>

namespace hemelb
{
  namespace topology
  {
    class Communicator
    {
      public:
        /**
         * Constructor for an uninitialised, NULL communicator
         * @param communicator
         */
        Communicator() :
            size(0)
        {

        }

        /**
         * Constructor to get data needed from an MPI communicator
         * @param communicator
         */
        Communicator(MPI_Comm communicator) :
            communicator(communicator)
        {
          int commRank, commSize;

          MPI_Comm_rank(communicator, &commRank);
          MPI_Comm_size(communicator, &commSize);
          MPI_Comm_group(communicator, &group);

          rank = commRank;
          size = commSize;
        }

        /**
         * Returns the local rank on the communicator
         * @return
         */
        inline proc_t GetRank() const
        {
          return rank;
        }

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        inline proc_t GetSize() const
        {
          return size;
        }

        /**
         * Returns the MPI communicator being used.
         * @return
         */
        inline MPI_Comm GetCommunicator() const
        {
          return communicator;
        }

        /**
         * Returns the MPI group being used.
         * @return
         */
        inline MPI_Group GetGroup() const
        {
          return group;
        }

      private:
        proc_t rank;
        proc_t size;
        MPI_Comm communicator;
        MPI_Group group;
    };
  }
}

#endif /* HEMELB_TOPOLOGY_COMMUNICATOR_H */
