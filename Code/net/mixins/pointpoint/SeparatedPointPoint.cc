// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "net/mixins/pointpoint/SeparatedPointPoint.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace net
  {

    void SeparatedPointPoint::EnsureEnoughRequests(size_t count)
    {
      if (requests.size() < count)
      {
        requests.resize(count, MPI_Request());
        statuses.resize(count, MPI_Status());
      }
    }

    void SeparatedPointPoint::ReceivePointToPoint()
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
                  &requests[m]);
        ++m;
      }

    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void SeparatedPointPoint::EnsurePreparedToSendReceive()
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

    void SeparatedPointPoint::SendPointToPoint()
    {

      // Make sure the MPI datatypes have been created.
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
                  &requests[receiveProcessorComms.size() + m]);

        ++m;
      }
    }

    /*!
     Free the allocated data.
     */
    SeparatedPointPoint::~SeparatedPointPoint()
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

    void SeparatedPointPoint::WaitPointToPoint()
    {

      MPI_Waitall((int) (sendProcessorComms.size() + receiveProcessorComms.size()), &requests[0], &statuses[0]);

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
      sendReceivePrepped = false;

    }
  }
}
