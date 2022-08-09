// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/pointpoint/CoalescePointPoint.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace net
  {

    void CoalescePointPoint::EnsureEnoughRequests(size_t count)
    {
      if (requests.size() < count)
      {
        requests.resize(count, MPI_Request());
        statuses.resize(count, MPI_Status());
      }
    }

    void CoalescePointPoint::ReceivePointToPoint()
    {

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (auto& [pid, pc]: receiveProcessorComms) {
        MPI_Irecv(pc.front().Pointer,
                  1,
                  pc.Type,
                  pid,
                  10,
                  communicator,
                  &requests[m]);
        ++m;
      }

      //if(m>1) {
      //  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("RecvPointToPoint() Neighbouring proc count: %i", m);
      //}

    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void CoalescePointPoint::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      for (auto& [_, pc]: sendProcessorComms)
      {
        pc.CreateMPIType();
      }

      for (auto& [_, pc]: receiveProcessorComms)
      {
        pc.CreateMPIType();
      }

      EnsureEnoughRequests(receiveProcessorComms.size() + sendProcessorComms.size());

      sendReceivePrepped = true;
    }

    void CoalescePointPoint::SendPointToPoint()
    {

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (auto& [pid, pc]: sendProcessorComms)
      {

        int TypeSizeStorage = 0; //DTMP:byte size tracking
        MPI_Type_size(pc.Type, &TypeSizeStorage); //DTMP:
        BytesSent += TypeSizeStorage; //DTMP:

        MPI_Isend(pc.front().Pointer,
                  1,
                  pc.Type,
                  pid,
                  10,
                  communicator,
                  &requests[receiveProcessorComms.size() + m]);

        ++m;
      }
    }

    /*!
     Free the allocated data.
     */
    CoalescePointPoint::~CoalescePointPoint()
    {
      if (sendReceivePrepped)
      {
        for (auto& [_, pc]: sendProcessorComms) {
          MPI_Type_free(&pc.Type);
        }

        for (auto& [_, pc]: receiveProcessorComms) {
          MPI_Type_free(&pc.Type);
        }

      }
    }

    void CoalescePointPoint::WaitPointToPoint()
    {

      MPI_Waitall((int) (sendProcessorComms.size() + receiveProcessorComms.size()),
                  requests.data(),
                  statuses.data());

      for (auto& [_, pc]: receiveProcessorComms) {
        MPI_Type_free(&pc.Type);
      }
      receiveProcessorComms.clear();

        for (auto& [_, pc]: sendProcessorComms) {
        MPI_Type_free(&pc.Type);
      }
      sendProcessorComms.clear();
      sendReceivePrepped = false;

    }
  }
}
