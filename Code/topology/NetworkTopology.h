#ifndef HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H
#define HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>

namespace hemelb
{
  namespace topology
  {
    class NeighbouringProcessor
    {
      public:
        // Rank of the neighbouring processor.
        unsigned int Rank;

        // The number of distributions shared between this neighbour and the current processor.
        int SharedFCount;

        // Index on this processor of the first distribution shared between this
        // neighbour and the current processor.
        int FirstSharedF;
    };

    class NetworkTopology
    {
      public:
        NetworkTopology(int * argCount, char *** argList, bool * oMachineDiscoverySuccess);

        ~NetworkTopology();

        bool IsCurrentProcTheIOProc() const;

        // Functions for getting the rank of this processor and the total size
        // of the topology.
        unsigned int GetLocalRank() const;
        unsigned int GetProcessorCount() const;
        int GetDepths() const;
        unsigned int GetMachineCount() const;

        // Number of local distributions shared with neighbouring processors.
        int TotalSharedFs;
        // The vector of all neighbouring processors.
        std::vector<NeighbouringProcessor> NeighbouringProcs;
        // For each processor in the topology, holds the index into the
        // neighbouring processor vector.
        short int * NeighbourIndexFromProcRank;
        // Array containing numbers of fluid sites on each processor.
        unsigned int * FluidSitesOnEachProcessor;

      private:
        bool InitialiseMachineInfo();

        unsigned int localRank;
        unsigned int processorCount;

        // Number of depths in the topology.
        int depths;
        // Number of machines in the topology.
        unsigned int machineCount;

        // Number of processors on each machine
        unsigned int * ProcCountOnEachMachine;
        // Machine Id where each processor is.
        unsigned int * MachineIdOfEachProc;
    };
  }
}

#endif /* HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H */
