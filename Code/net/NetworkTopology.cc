// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/NetworkTopology.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace net
  {
    NetworkTopology NetworkTopology::instance;
    bool NetworkTopology::initialised = false;

    // Must be specified to prevent the default constructor being public.
    NetworkTopology::NetworkTopology()
    {

    }

    NetworkTopology* NetworkTopology::Instance()
    {
      return &instance;
    }

    void NetworkTopology::Init(MpiCommunicator& commun)
    {
      if (!initialised)
      {
        initialised = true;
        comms = commun;
        InitialiseMachineInfo();
      }
    }

    NetworkTopology::~NetworkTopology()
    {
      if (initialised)
      {
        delete[] ProcCountOnEachMachine;
        delete[] MachineIdOfEachProc;
      }
    }

    bool NetworkTopology::IsCurrentProcTheIOProc() const
    {
      return GetLocalRank() == GetIOProcRank();
    }

    int NetworkTopology::GetIOProcRank() const
    {
      return 0;
    }

    proc_t NetworkTopology::GetLocalRank() const
    {
      return comms.Rank();
    }

    proc_t NetworkTopology::GetProcessorCount() const
    {
      return comms.Size();
    }

    int NetworkTopology::GetDepths() const
    {
      return depths;
    }

    unsigned int NetworkTopology::GetMachineCount() const
    {
      return machineCount;
    }

  }
}
