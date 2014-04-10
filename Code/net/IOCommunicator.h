// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_IOCOMMUNICATOR_H
#define HEMELB_NET_IOCOMMUNICATOR_H

//#include <vector>
//#include <cstdio>
//
//#include "constants.h"
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {
    /**
     * An MPI communicator with a special I/O rank.
     */
    class IOCommunicator : public MpiCommunicator
    {
      public:
        IOCommunicator(const MpiCommunicator& comm);
        bool OnIORank() const;
        int GetIORank() const;
    };
  }
}

#endif /* HEMELB_NET_IOCOMMUNICATOR_H */
