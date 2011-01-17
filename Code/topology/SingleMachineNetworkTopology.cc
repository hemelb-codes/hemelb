#include "topology/NetworkTopology.h"

/*
 * NOTE!
 *
 * Changes to this file may need to be made in concert with changes to
 * MultiMachineNetworkTopology.cc
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
      // the machine is assumed to be only one if this function is
      // used instead of the previous one

      depths = 1;
      machineCount = 1;

      MachineIdOfEachProc = new int[GetProcessorCount()];
      ProcCountOnEachMachine = new int[machineCount];

      for (int i = 0; i < GetProcessorCount(); i++)
      {
        MachineIdOfEachProc[i] = 0;
      }
      ProcCountOnEachMachine[0] = GetProcessorCount();

      return true;
    }
  }
}

