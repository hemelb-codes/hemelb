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

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); ++it)
      {
        MPI_Irecv(it->second.PointerList.front(),
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  communicator,
                  &mRequests[m]);
        ++m;
      }
    }

    void Net::Send()
    {
      // Make sure the datatypes have been created.
      EnsurePreparedToSendReceive();

      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); ++it)
      {
        int TypeSizeStorage = 0;                         //DTMP:byte size tracking 
        MPI_Type_size(it->second.Type,&TypeSizeStorage); //DTMP:
        BytesSent += TypeSizeStorage;                   //DTMP:


        MPI_Isend(it->second.PointerList.front(),
                  1,
                  it->second.Type,
                  it->first,
                  10,
                  communicator,
                  &mRequests[mReceiveProcessorComms.size() + m]);

        ++m;
      }
    }

    void Net::Wait()
    {
      SyncPointsCounted++; //DTMP: counter for monitoring purposes.

      MPI_Waitall((int) (mSendProcessorComms.size() + mReceiveProcessorComms.size()),
                  &mRequests[0],
                  &mStatuses[0]);

      sendReceivePrepped = false;

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      mReceiveProcessorComms.clear();

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      mSendProcessorComms.clear();
    }

    // Helper function to get the ProcessorCommunications object, and create it if it doesn't exist yet.
    Net::ProcComms* Net::GetProcComms(proc_t iRank, bool iIsSend)
    {
      std::map<proc_t, ProcComms>* lMap = iIsSend
        ? &mSendProcessorComms
        : &mReceiveProcessorComms;

      std::map<proc_t, ProcComms>::iterator lValue = lMap->find(iRank);

      if (lValue == lMap->end())
      {
        ProcComms lRet;
        lMap->insert(std::pair<proc_t, ProcComms>(iRank, lRet));
        return &lMap->find(iRank)->second;
      }
      else
      {
        return & (lValue ->second);
      }
    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void Net::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); ++it)
      {
        CreateMPIType(& (*it).second);
      }

      for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); ++it)
      {
        CreateMPIType(& (*it).second);
      }

      EnsureEnoughRequests(mReceiveProcessorComms.size() + mSendProcessorComms.size());

      sendReceivePrepped = true;
    }

    // Helper function to create a MPI derived datatype given a list of pointers, types and lengths.
    void Net::CreateMPIType(ProcComms *iMetaData)
    {
      MPI_Aint* displacements = new MPI_Aint[iMetaData->PointerList.size()];
      int* lengths = new int[iMetaData->PointerList.size()];
      MPI_Datatype* types = new MPI_Datatype[iMetaData->PointerList.size()];

      int lLocation = 0;

      for (std::vector<void*>::const_iterator it = iMetaData->PointerList.begin(); it
          != iMetaData->PointerList.end(); ++it)
      {
        MPI_Get_address(*it, &displacements[lLocation]);
        ++lLocation;
      }

      for (int ii = (int) iMetaData->PointerList.size() - 1; ii >= 0; ii--)
      {
        displacements[ii] -= displacements[0];

        lengths[ii] = iMetaData->LengthList[ii];
        types[ii] = iMetaData->TypeList[ii];
        /* total_length == lengths[ii]*/
      }

      // Create the type and commit it.
      MPI_Type_create_struct((int) iMetaData->PointerList.size(),
                             lengths,
                             displacements,
                             types,
                             &iMetaData->Type);
      MPI_Type_commit(&iMetaData->Type);

      delete[] displacements;
      delete[] lengths;
      delete[] types;
    }

    Net::Net(): 
        BytesSent(0), 
        SyncPointsCounted(0)
    {
      sendReceivePrepped = false;
      communicator = MPI_COMM_WORLD;
    }

    Net::Net(MPI_Comm commObject)
    {
      sendReceivePrepped = false;
      communicator = commObject;
    }

    /*!
     Free the allocated data.
     */
    Net::~Net()
    {
      if (sendReceivePrepped)
      {
        for (std::map<proc_t, ProcComms>::iterator it = mSendProcessorComms.begin(); it
            != mSendProcessorComms.end(); ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

        for (std::map<proc_t, ProcComms>::iterator it = mReceiveProcessorComms.begin(); it
            != mReceiveProcessorComms.end(); ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

      }
    }

  }
}
