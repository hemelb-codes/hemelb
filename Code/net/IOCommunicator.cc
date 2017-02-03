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
    IOCommunicator::IOCommunicator(const MpiCommunicator& comm) :
        MpiCommunicator(comm)
    {
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
