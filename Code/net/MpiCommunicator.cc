//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace net
  {

    MpiCommunicator::MpiCommunicator() :
      size(0)
    {

    }

    MpiCommunicator::MpiCommunicator(MPI_Comm communicator) :
      communicator(communicator)
    {
      int commRank, commSize;

      MPI_Comm_rank(communicator, &commRank);
      MPI_Comm_size(communicator, &commSize);
      MPI_Comm_group(communicator, &group);

      rank = commRank;
      size = commSize;
    }

    MpiCommunicator::MpiCommunicator(proc_t rank, proc_t size) :
      rank(rank), size(size), communicator(NULL), group(NULL)
    {
    }

  }
}
