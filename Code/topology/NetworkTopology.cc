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
      MPI_Comm_size(MPI_COMM_WORLD, &processorCount);
      MPI_Comm_rank(MPI_COMM_WORLD, &localRank);

      if (IsCurrentProcTheIOProc())
      {
        printf("thread_level_provided %i\n", thread_level_provided);
      }
    }

    NetworkTopology::~NetworkTopology()
    {
      MPI_Finalize();

      delete[] NeighbourIndexFromProcRank;
      delete[] FluidSitesOnEachProcessor;
      delete[] ProcCountOnEachMachine;
      delete[] MachineIdOfEachProc;
    }

    bool NetworkTopology::IsCurrentProcTheIOProc() const
    {
      return localRank == 0;
    }

    int NetworkTopology::GetLocalRank() const
    {
      return localRank;
    }

    int NetworkTopology::GetProcessorCount() const
    {
      return processorCount;
    }

  }
}
