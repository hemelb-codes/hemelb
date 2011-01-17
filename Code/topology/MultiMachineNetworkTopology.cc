#include "topology/NetworkTopology.h"

/*
 * NOTE!
 *
 * Changes to this file may need to be made in concert with changes to
 * SingleMachineNetworkTopology.cc
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
    bool NetworkTopology::InitialiseMachineInfo()
    {
      int err;
      int *depth, **color;
      int machine_id, flag, is_found;
      int i, j, sum;

      Depths = 0;

      err = MPI_Attr_get(MPI_COMM_WORLD, MPICHX_TOPOLOGY_DEPTHS, &depth, &flag);

      if (err != MPI_SUCCESS || flag == 0)
      {
        return false;
      }

      err = MPI_Attr_get(MPI_COMM_WORLD, MPICHX_TOPOLOGY_COLORS, &color, &flag);

      if (err != MPI_SUCCESS || flag == 0)
      {
        return false;
      }

      MachineCount = 0;

      MachineIdOfEachProc = new int[GetProcessorCount()];
      ProcCountOnEachMachine = new int[GetProcessorCount()];

      for (i = 0; i < GetProcessorCount(); i++)
      {
        ProcCountOnEachMachine[i] = 0;
      }
      for (i = 0; i < GetProcessorCount(); i++)
      {
        if (depth[i] != 4)
          continue;

        Depths = max(Depths, depth[i]);

        for (j = 0, is_found = 0; j < MachineCount && is_found == 0; j++)
        {
          if (color[i][3] == MachineIdOfEachProc[j])
          {
            is_found = 1;
            ++ProcCountOnEachMachine[MachineIdOfEachProc[j]];
          }
        }
        if (is_found == 1)
          continue;

        MachineIdOfEachProc[MachineCount] = color[i][3];
        ++ProcCountOnEachMachine[MachineCount];
        ++MachineCount;
      }
      MachineCount = max(1, MachineCount);

      if (MachineCount == 1)
      {
        for (i = 0; i < GetProcessorCount(); i++)
        {
          MachineIdOfEachProc[i] = 0;
        }
        ProcCountOnEachMachine[0] = GetProcessorCount();
      }
      else
      {
        for (i = 0; i < GetProcessorCount(); i++)
        {
          sum = 0;
          machine_id = 0;

          is_found = 0;

          while (!is_found)
          {
            if (sum + ProcCountOnEachMachine[machine_id] > i)
            {
              is_found = 1;
              continue;
            }
            sum += ProcCountOnEachMachine[machine_id];
            ++machine_id;
          }
          MachineIdOfEachProc[i] = machine_id;
        }
      }
      return true;
    }
  }
}
