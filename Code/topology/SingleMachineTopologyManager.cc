#include "topology/TopologyManager.h"

/*
 * NOTE!
 *
 * Changes to this file may need to be made in concert with changes to
 * MultiMachineTopologyManager.cc
 *
 * Only one of these two files will be compiled on any given build, and
 * they determine whether our topology will be over multiple machines or
 * not.
 */

namespace hemelb
{
  namespace topology
  {

    /*!
     If one has more than one machine. The topology discovery mechanism is implemented in this function
     */
    TopologyManager::TopologyManager(NetworkTopology* bNetworkTopology,
                                     bool* oSuccess)
    {
      // the machine is assumed to be only one if this function is
      // used instead of the previous one

      bNetworkTopology->Depths = 1;
      bNetworkTopology->MachineCount = 1;

      bNetworkTopology->MachineIdOfEachProc
          = new int[bNetworkTopology->GetProcessorCount()];
      bNetworkTopology->ProcCountOnEachMachine
          = new int[bNetworkTopology->MachineCount];

      for (int i = 0; i < bNetworkTopology->GetProcessorCount(); i++)
      {
        bNetworkTopology->MachineIdOfEachProc[i] = 0;
      }
      bNetworkTopology->ProcCountOnEachMachine[0]
          = bNetworkTopology->GetProcessorCount();

      *oSuccess = true;
    }
  }
}

