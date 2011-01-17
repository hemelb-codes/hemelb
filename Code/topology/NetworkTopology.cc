#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace topology
  {

    NeighbouringProcessor::~NeighbouringProcessor()
    {
      delete[] SharedFReceivingIndex;
    }

    NetworkTopology::NetworkTopology(int * argCount, char *** argList)
    {
      int thread_level_provided;

      MPI_Init_thread(argCount, argList, MPI_THREAD_FUNNELED,
                      &thread_level_provided);
      MPI_Comm_size(MPI_COMM_WORLD, &ProcessorCount);
      MPI_Comm_rank(MPI_COMM_WORLD, &LocalRank);

      if (IsCurrentProcTheIOProc())
      {
        printf("thread_level_provided %i\n", thread_level_provided);
      }
    }

    NetworkTopology::~NetworkTopology()
    {
      delete[] NeighbourIndexFromProcRank;
      delete[] FluidSitesOnEachProcessor;
      delete[] ProcCountOnEachMachine;
      delete[] MachineIdOfEachProc;
    }

    bool NetworkTopology::IsCurrentProcTheIOProc() const
    {
      return LocalRank == 0;
    }

  }
}
