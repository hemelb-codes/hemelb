#ifndef HEMELB_TOPOLOGY_NETWORKTOPOLOGY_H
#define HEMELB_TOPOLOGY_NETWORKTOPOLOGY_H

#include <vector>
#include <cstdio>

#include "constants.h"
#include "topology/Communicator.h"

namespace hemelb
{
  namespace topology
  {
    class NetworkTopology
    {
      public:
        static NetworkTopology* Instance();
        ~NetworkTopology();

        void Init(int argCount, char ** argList, bool * oMachineDiscoverySuccess);
        bool IsCurrentProcTheIOProc() const;

        // Functions for getting the rank of this processor and the total size
        // of the topology.
        proc_t GetLocalRank() const;
        proc_t GetProcessorCount() const;
        int GetDepths() const;
        unsigned int GetMachineCount() const;

      private:
        NetworkTopology();
        bool InitialiseMachineInfo();

        Communicator comms;

        // Number of depths in the topology.
        int depths;
        // Number of machines in the topology.
        unsigned int machineCount;

        // Number of processors on each machine
        proc_t* ProcCountOnEachMachine;
        // Machine Id where each processor is.
        unsigned int* MachineIdOfEachProc;

        /**
         * This variable is necessary, because the destructor for this static object will always
         * be called, regardless of whether the init method (that actually initialises the MPI
         * environment) is called.
         */
        static bool initialised;
        static NetworkTopology instance;
    };
  }
}

#endif /* HEMELB_TOPOLOGY_NETWORKTOPOLOGY_H */
