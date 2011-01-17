#ifndef HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H
#define HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>
#include "mpiInclude.h"
#include "lb/GlobalLatticeData.h"

namespace hemelb
{
  namespace topology
  {
    class NeighbouringProcessor
    {
      public:
        ~NeighbouringProcessor();

        // Rank of the neighbouring processor.
        int Rank;

        // The number of distributions shared between this neighbour and the current processor.
        int SharedFCount;

        // Index on this processor of the first distribution shared between this
        // neighbour and the current processor.
        int FirstSharedF;

        // Array of the indexes corresponding to received distributions.
        int * SharedFReceivingIndex;
    };

    class NetworkTopology
    {
      public:
        NetworkTopology(int * argCount,
                        char *** argList,
                        bool * oMachineDiscoverySuccess);

        ~NetworkTopology();

        bool IsCurrentProcTheIOProc() const;

        // Functions for getting the rank of this processor and the total size
        /// of the topology.
        int GetLocalRank() const;
        int GetProcessorCount() const;

        void DecomposeDomain(int iTotalFluidSites,
                             const lb::GlobalLatticeData & bGlobLatDat);

        // Number of local distributions shared with neighbouring processors.
        int TotalSharedFs;
        // The vector of all neighbouring processors.
        std::vector<NeighbouringProcessor*> NeighbouringProcs;
        // For each processor in the topology, holds the index into the
        // neighbouring processor vector.
        short int * NeighbourIndexFromProcRank;
        // Array containing numbers of fluid sites on each processor.
        int * FluidSitesOnEachProcessor;

        // Number of processors on each machine
        int * ProcCountOnEachMachine;
        // Machine Id where each processor is.
        int * MachineIdOfEachProc;
        // Number of depths in the topology.
        int Depths;
        // Number of machines in the topology.
        int MachineCount;

      private:
        struct SiteLocation
        {
            short int i, j, k;
        };

        bool InitialiseMachineInfo();

        void
        AssignFluidSitesToProcessors(int & proc_count,
                                     int & fluid_sites_per_unit,
                                     int & unvisited_fluid_sites,
                                     const int iCurrentProcId,
                                     const bool iIsMachineLevel,
                                     const lb::GlobalLatticeData &iGlobLatDat);

        int localRank;
        int processorCount;
    };
  }
}

#endif /* HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H */
