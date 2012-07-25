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
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        for (ProcComms::iterator request = it->second.begin(); request != it->second.end(); request++)
        {

          MPI_Irecv(request->Pointer,
                    request->Count,
                    request->Type,
                    it->first,
                    10,
                    communicator.GetCommunicator(),
                    &requests[m]);
          ++m;
        }
      }

    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void SeparatedPointPoint::EnsurePreparedToSendReceive()
    {
      count_sends=0;
      count_receives=0;
      if (sendReceivePrepped)
      {
        return;
      }

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        count_sends += it->second.size();
      }

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        count_receives += it->second.size();
      }

      EnsureEnoughRequests(count_sends+count_receives);

      sendReceivePrepped = true;
    }

    void SeparatedPointPoint::SendPointToPoint()
    {

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        for (ProcComms::iterator request = it->second.begin(); request != it->second.end(); request++)
        {
          MPI_Isend(request->Pointer,
                    request->Count,
                    request->Type,
                    it->first,
                    10,
                    communicator.GetCommunicator(),
                    &requests[count_receives+m]);
          ++m;
        }
      }
    }

    /*!
     Free the allocated data.
     */
    SeparatedPointPoint::~SeparatedPointPoint()
    {
    }

    void SeparatedPointPoint::WaitPointToPoint()
    {
      MPI_Waitall(static_cast<int>(count_sends+count_receives), &requests[0], &statuses[0]);

      receiveProcessorComms.clear();
      sendProcessorComms.clear();
      sendReceivePrepped = false;

    }
  }
}
