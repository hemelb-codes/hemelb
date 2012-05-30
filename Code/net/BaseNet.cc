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

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();

      ReceiveGathers();
      ReceiveGatherVs();
      ReceivePointToPoint();
    }

    void BaseNet::Send()
    {
      // Make sure the datatypes have been created.
      EnsurePreparedToSendReceive();
      SendPointToPoint();
      SendGathers();
      SendGatherVs();

    }

    void BaseNet::Wait()
    {
      SyncPointsCounted++; //DTMP: counter for monitoring purposes.
      WaitPointToPoint();
      WaitGathers();
      WaitGatherVs();
    }

  }
}
