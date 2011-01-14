#include "topology/TopologyManager.h"

/*
 * NOTE!
 *
 * Changes to this file may need to be made in concert with changes to
 * SingleMachineTopologyManager.cc
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
    void TopologyManager::FindTopology(NetworkTopology* bNetworkTopology,
                                       bool & oWasSuccessful)
    {
      int err;
      int *depth, **color;
      int machine_id, flag, is_found;
      int i, j, sum;

      bNetworkTopology->Depths = 0;

      err = MPI_Attr_get(MPI_COMM_WORLD, MPICHX_TOPOLOGY_DEPTHS, &depth, &flag);

      if (err != MPI_SUCCESS || flag == 0)
      {
        oWasSuccessful = false;
        return;
      }

      err = MPI_Attr_get(MPI_COMM_WORLD, MPICHX_TOPOLOGY_COLORS, &color, &flag);

      if (err != MPI_SUCCESS || flag == 0)
      {
        oWasSuccessful = false;
        return;
      }

      bNetworkTopology->MachineCount = 0;

      bNetworkTopology->MachineIdOfEachProc
          = new int[bNetworkTopology->ProcessorCount];
      bNetworkTopology->ProcCountOnEachMachine
          = new int[bNetworkTopology->ProcessorCount];

      for (i = 0; i < bNetworkTopology->ProcessorCount; i++)
      {
        bNetworkTopology->ProcCountOnEachMachine[i] = 0;
      }
      for (i = 0; i < bNetworkTopology->ProcessorCount; i++)
      {
        if (depth[i] != 4)
          continue;

        bNetworkTopology->Depths = max(bNetworkTopology->Depths, depth[i]);

        for (j = 0, is_found = 0; j < bNetworkTopology->MachineCount
            && is_found == 0; j++)
        {
          if (color[i][3] == bNetworkTopology->MachineIdOfEachProc[j])
          {
            is_found = 1;
            ++bNetworkTopology->ProcCountOnEachMachine[bNetworkTopology->MachineIdOfEachProc[j]];
          }
        }
        if (is_found == 1)
          continue;

        bNetworkTopology->MachineIdOfEachProc[bNetworkTopology->MachineCount]
            = color[i][3];
        ++bNetworkTopology->ProcCountOnEachMachine[bNetworkTopology->MachineCount];
        ++bNetworkTopology->MachineCount;
      }
      bNetworkTopology->MachineCount = max(1, bNetworkTopology->MachineCount);

      if (bNetworkTopology->MachineCount == 1)
      {
        for (i = 0; i < bNetworkTopology->ProcessorCount; i++)
        {
          bNetworkTopology->MachineIdOfEachProc[i] = 0;
        }
        bNetworkTopology->ProcCountOnEachMachine[0]
            = bNetworkTopology->ProcessorCount;
      }
      else
      {
        for (i = 0; i < bNetworkTopology->ProcessorCount; i++)
        {
          sum = 0;
          machine_id = 0;

          is_found = 0;

          while (!is_found)
          {
            if (sum + bNetworkTopology->ProcCountOnEachMachine[machine_id] > i)
            {
              is_found = 1;
              continue;
            }
            sum += bNetworkTopology->ProcCountOnEachMachine[machine_id];
            ++machine_id;
          }
          bNetworkTopology->MachineIdOfEachProc[i] = machine_id;
        }
      }
      oWasSuccessful = true;
    }
  }
}
