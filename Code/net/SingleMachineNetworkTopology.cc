// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/NetworkTopology.h"

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
  namespace net
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

      MachineIdOfEachProc = new unsigned int[GetProcessorCount()];
      ProcCountOnEachMachine = new proc_t[machineCount];

      for (proc_t i = 0; i < GetProcessorCount(); i++)
      {
        MachineIdOfEachProc[i] = 0;
      }
      ProcCountOnEachMachine[0] = GetProcessorCount();

      return true;
    }
  }
}

