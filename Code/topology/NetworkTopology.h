#ifndef HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H
#define HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>

#include "constants.h"

namespace hemelb
{
  namespace topology
  {
    class NeighbouringProcessor
    {
      public:
        // Rank of the neighbouring processor.
        proc_t Rank;

        // The number of distributions shared between this neighbour and the current processor.
        site_t SharedFCount;

        // Index on this processor of the first distribution shared between this
        // neighbour and the current processor.
        site_t FirstSharedF;
    };

    class NetworkTopology
    {
      public:
        static NetworkTopology* Instance();

        void Init(int * argCount, char *** argList, bool * oMachineDiscoverySuccess);
        bool IsCurrentProcTheIOProc() const;

        // Functions for getting the rank of this processor and the total size
        // of the topology.
        proc_t GetLocalRank() const;
        proc_t GetProcessorCount() const;
        int GetDepths() const;
        unsigned int GetMachineCount() const;

        // Number of local distributions shared with neighbouring processors.
        site_t TotalSharedFs;
        // The vector of all neighbouring processors.
        std::vector<NeighbouringProcessor> NeighbouringProcs;
        // For each processor in the topology, holds the index into the
        // neighbouring processor vector.
        proc_t* NeighbourIndexFromProcRank;
        // Array containing numbers of fluid sites on each processor.
        site_t* FluidSitesOnEachProcessor;

      private:
        NetworkTopology();
        ~NetworkTopology();
        bool InitialiseMachineInfo();

        proc_t localRank;
        proc_t processorCount;

        // Number of depths in the topology.
        int depths;
        // Number of machines in the topology.
        unsigned int machineCount;

        // Number of processors on each machine
        proc_t* ProcCountOnEachMachine;
        // Machine Id where each processor is.
        unsigned int* MachineIdOfEachProc;

        static bool initialised;
        static NetworkTopology instance;
    };
  }
}

#endif /* HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H */
