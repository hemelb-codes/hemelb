// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

      for (auto& [pid, pc]: receiveProcessorComms)
      {
        for (auto& req: pc)
        {

          MPI_Irecv(req.Pointer,
                    req.Count,
                    req.Type,
                    pid,
                    10,
                    communicator,
                    &requests[m]);
          ++m;
        }
      }

    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void SeparatedPointPoint::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      count_sends = 0;
      count_receives = 0;
      for (auto& [_, pc]: sendProcessorComms)
      {
        count_sends += pc.size();
      }

      for (auto& [_, pc]: receiveProcessorComms)
      {
        count_receives += pc.size();
      }

      EnsureEnoughRequests(count_sends + count_receives);

      sendReceivePrepped = true;
    }

    void SeparatedPointPoint::SendPointToPoint()
    {

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (auto& [pid, pc]: sendProcessorComms)
      {
        for (auto& req: pc)
        {
          MPI_Isend(req.Pointer,
                    req.Count,
                    req.Type,
                    pid,
                    10,
                    communicator,
                    &requests[count_receives + m]);
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
      MPI_Waitall(static_cast<int>(count_sends + count_receives), &requests[0], &statuses[0]);

      receiveProcessorComms.clear();
      sendProcessorComms.clear();
      sendReceivePrepped = false;

    }
  }
}
