/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "net/net.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace net
  {

    void Net::EnsureEnoughRequests(size_t count)
    {
      if (mRequests.size() < count)
      {
        size_t deficit = count - mRequests.size();
        for (unsigned int ii = 0; ii < deficit; ii++)
        {
          mRequests.push_back(MPI_Request());
          mStatuses.push_back(MPI_Status());
        }
      }
    }

    void Net::Dispatch()
    {
      Send();
      Receive();
      Wait();
    }

    /*!
     This is called from the main function.  First function to deal with processors.
     The domain partitioning technique and the management of the
     buffers useful for the inter-processor communications are
     implemented in this function.  The domain decomposition is based
     on a graph growing partitioning technique.
     */

    void Net::Receive()
    {
      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();

      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        MPI_Irecv(it->second.front().Pointer,
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  communicator.GetCommunicator(),
                  &mRequests[m]);
        ++m;
      }

      int gather_index = 0;
      for (std::vector<ScalarRequest>::iterator it = gatherReceiveProcessorComms.begin();
          it != gatherReceiveProcessorComms.end(); ++it)
      {
        ScalarRequest toself = gatherSendProcessorComms[communicator.GetRank()][m];
        MPI_Gather(toself.Pointer,
                   1,
                   toself.Type,
                   it->Pointer,
                   1,
                   it->Type,
                   communicator.GetRank(),
                   communicator.GetCommunicator());
        ++gather_index;
      }
    }

    void Net::Send()
    {
      // Make sure the datatypes have been created.
      EnsurePreparedToSendReceive();

      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        int TypeSizeStorage = 0; //DTMP:byte size tracking
        MPI_Type_size(it->second.Type, &TypeSizeStorage); //DTMP:
        BytesSent += TypeSizeStorage; //DTMP:

        MPI_Isend(it->second.front().Pointer,
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  communicator.GetCommunicator(),
                  &mRequests[receiveProcessorComms.size() + m]);

        ++m;
      }

      for (std::map<proc_t, GatherProcComms>::iterator it = gatherSendProcessorComms.begin();
          it != gatherSendProcessorComms.end(); ++it)
      {
        if (it->first==communicator.GetRank())
        {
          continue;
        }
        for (std::vector<ScalarRequest>::iterator req = it->second.begin(); req != it->second.end(); req++)
        {
          MPI_Gather(req->Pointer,
                     1,
                     req->Type,
                     NULL,
                     1,
                     req->Type,
                     it->first,
                     communicator.GetCommunicator());
        }

      }
    }

    void Net::Wait()
    {
      SyncPointsCounted++; //DTMP: counter for monitoring purposes.

      MPI_Waitall((int) (sendProcessorComms.size() + receiveProcessorComms.size()), &mRequests[0], &mStatuses[0]);

      sendReceivePrepped = false;

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      receiveProcessorComms.clear();

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      sendProcessorComms.clear();
    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void Net::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        it->second.CreateMPIType();
      }

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        it->second.CreateMPIType();
      }

      EnsureEnoughRequests(receiveProcessorComms.size() + sendProcessorComms.size());

      sendReceivePrepped = true;
    }

    Net::Net() :
        BytesSent(0), SyncPointsCounted(0), communicator(topology::NetworkTopology::Instance()->GetComms())
    {
      sendReceivePrepped = false;
    }

    Net::Net(topology::Communicator &commObject) :
        BytesSent(0), SyncPointsCounted(0), communicator(commObject)
    {
      sendReceivePrepped = false;
    }

    /*!
     Free the allocated data.
     */
    Net::~Net()
    {
      if (sendReceivePrepped)
      {
        for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end();
            ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

        for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
            it != receiveProcessorComms.end(); ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

      }
    }

  }
}
