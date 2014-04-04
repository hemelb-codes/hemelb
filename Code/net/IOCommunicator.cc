// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/IOCommunicator.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace net
  {
    IOCommunicator IOCommunicator::instance;
    bool IOCommunicator::initialised = false;

    // Must be specified to prevent the default constructor being public.
    IOCommunicator::IOCommunicator()
    {

    }

    IOCommunicator* IOCommunicator::Instance()
    {
      return &instance;
    }

    void IOCommunicator::Init(MpiCommunicator& commun)
    {
      if (!initialised)
      {
        initialised = true;
        comms = commun;
      }
    }

    IOCommunicator::~IOCommunicator()
    {
    }

    bool IOCommunicator::IsCurrentProcTheIOProc() const
    {
      return Rank() == GetIOProcRank();
    }

    int IOCommunicator::GetIOProcRank() const
    {
      return 0;
    }

    proc_t IOCommunicator::Rank() const
    {
      return comms.Rank();
    }

    proc_t IOCommunicator::Size() const
    {
      return comms.Size();
    }

  }
}
