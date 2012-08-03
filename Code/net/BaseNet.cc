// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "net/BaseNet.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace net
  {

    void BaseNet::Dispatch()
    {
      Send();
      Receive();
      Wait();
    }

    BaseNet::BaseNet() :
        BytesSent(0), SyncPointsCounted(0), communicator(topology::NetworkTopology::Instance()->GetComms())
    {
    }

    BaseNet::BaseNet(topology::Communicator &commObject) :
        BytesSent(0), SyncPointsCounted(0), communicator(commObject)
    {
    }

    void BaseNet::Receive()
    {
      ReceiveGathers();
      ReceiveGatherVs();
      ReceiveAllToAll();
      // Ensure collectives are called before point-to-point, as some implementing mixins implement collectives via point-to-point
      ReceivePointToPoint();
    }

    void BaseNet::Send()
    {

      SendGathers();
      SendGatherVs();
      SendAllToAll();
      // Ensure collectives are called before point-to-point, as some implementing mixins implement collectives via point-to-point
      SendPointToPoint();
    }

    void BaseNet::Wait()
    {
      SyncPointsCounted++; //DTMP: counter for monitoring purposes.

      WaitGathers();
      WaitGatherVs();
      WaitPointToPoint();
      WaitAllToAll();

      displacementsBuffer.clear();
      countsBuffer.clear();
    }

    std::vector<int> & BaseNet::GetDisplacementsBuffer()
    {
      displacementsBuffer.push_back(std::vector<int>());
      return displacementsBuffer.back();
    }

    std::vector<int> & BaseNet::GetCountsBuffer()
    {
      countsBuffer.push_back(std::vector<int>());
      return countsBuffer.back();
    }

  }
}
