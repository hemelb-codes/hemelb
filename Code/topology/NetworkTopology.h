#ifndef HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H
#define HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>

namespace hemelb
{
  namespace topology
  {
    struct NeighbouringProcessor
    {
        ~NeighbouringProcessor()
        {
          delete[] SharedFLocation;
          delete[] SharedFReceivingIndex;
        }

        // Rank of the neighbouring processor.
        int Rank;

        // The number of distributions shared between this neighbour and the current processor.
        int SharedFCount;

        // Site co-ordinates for each of the shared distribution, and the number of
        // the corresponding direction vector. Array is 4 elements for each shared distribution.
        short int * SharedFLocation;

        // Index on this processor of the first distribution shared between this
        // neighbour and the current processor.
        int FirstSharedF;

        // Array of the indexes corresponding to received distributions.
        int * SharedFReceivingIndex;
    };

    struct NetworkTopology
    {
        ~NetworkTopology()
        {
          delete[] NeighbourIndexFromProcRank;
          delete[] FluidSitesOnEachProcessor;
          delete[] ProcCountOnEachMachine;
          delete[] MachineIdOfEachProc;
        }

        bool IsCurrentProcTheIOProc() const
        {
          return LocalRank == 0;
        }

        char * GetCurrentProcIdentifier() const
        {
          // Note that this has enough space for an integer of 9 digits, so we're good up to a billion
          // processors.
          char* lRet = new char[15];
          sprintf(lRet, "rank %i", LocalRank);
          return lRet;
        }

        // The rank of this processor, the total size of the topology, and the
        // number of fluid sites on this rank.
        int LocalRank;
        int ProcessorCount;

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
    };
  }
}

#endif /* HEMELB_TOPOLOGY_LOCALNETWORKTOPOLOGY_H */
