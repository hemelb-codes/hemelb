// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/IOCommunicator.h"
#include "net/mpi.h"

namespace hemelb::net
{
    IOCommunicator::IOCommunicator(const MpiCommunicator& comm) :
        MpiCommunicator(comm),
        nodeComm(comm.SplitType()),
        amNodeLeader(nodeComm.Rank() == IO_RANK),
        leadersComm(comm.Split(amNodeLeader))
    {
    }

}
