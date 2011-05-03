#include "topology/NetworkTopology.h"
#include "math.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace topology
  {
    NetworkTopology::NetworkTopology(int * argCount, char *** argList, bool* oSuccess)
    {
      MPI_Init(argCount, argList);

      int tempSize = 0, tempRank = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &tempSize);
      MPI_Comm_rank(MPI_COMM_WORLD, &tempRank);

      processorCount = (proc_t) tempSize;
      localRank = (proc_t) tempRank;

      *oSuccess = InitialiseMachineInfo();

      FluidSitesOnEachProcessor = new site_t[processorCount];
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

    proc_t NetworkTopology::GetLocalRank() const
    {
      return localRank;
    }

    proc_t NetworkTopology::GetProcessorCount() const
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
