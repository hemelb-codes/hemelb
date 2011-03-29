#include "topology/NetworkTopology.h"
#include "math.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace topology
  {
    NetworkTopology::NetworkTopology(int * argCount, char *** argList, bool* oSuccess)
    {
      int thread_level_provided = 0;

      MPI_Init_thread(argCount, argList, MPI_THREAD_FUNNELED, &thread_level_provided);

      int tempSize = 0, tempRank = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &tempSize);
      MPI_Comm_rank(MPI_COMM_WORLD, &tempRank);

      processorCount = tempSize;
      localRank = tempRank;

      if (IsCurrentProcTheIOProc())
      {
        printf("thread_level_provided %i\n", thread_level_provided);
      }

      *oSuccess = InitialiseMachineInfo();

      FluidSitesOnEachProcessor = new unsigned int[processorCount];
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

    unsigned int NetworkTopology::GetLocalRank() const
    {
      return localRank;
    }

    unsigned int NetworkTopology::GetProcessorCount() const
    {
      return processorCount;
    }

    int NetworkTopology::GetDepths() const
    {
      return depths;
    }

    unsigned int NetworkTopology::GetMachineCount() const
    {
      return machineCount;
    }

  }
}
