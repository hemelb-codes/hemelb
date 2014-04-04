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
    IOCommunicator::IOCommunicator() : MpiCommunicator()
    {

    }
    IOCommunicator::IOCommunicator(const MpiCommunicator& comm) : MpiCommunicator(comm)
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
        instance = commun;
      }
    }


    bool IOCommunicator::OnIORank() const
    {
      return Rank() == GetIORank();
    }

    int IOCommunicator::GetIORank() const
    {
      return 0;
    }

  }
}
