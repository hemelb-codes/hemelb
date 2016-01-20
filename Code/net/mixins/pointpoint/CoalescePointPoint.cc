
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

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {

        MPI_Irecv(it->second.front().Pointer,
                  1,
                  it->second.Type,
                  it->first,
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

    void CoalescePointPoint::SendPointToPoint()
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

    void CoalescePointPoint::WaitPointToPoint()
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
