#include "net/mixins/CoalescePointPoint.h"
#include "log/Logger.h"
namespace hemelb{
  namespace net{
    void CoalescePointPoint::EnsureEnoughRequests(size_t count)
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

    void CoalescePointPoint::ReceivePointToPoint()
    {

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
    }

    /*!
     Free the allocated data.
     */
    CoalescePointPoint::~CoalescePointPoint()
    {
      if (sendReceivePrepped)
      {
        for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

        for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
            ++it)
        {
          MPI_Type_free(&it->second.Type);
        }

      }
    }

    void CoalescePointPoint::Wait()
    {
      BaseNet::Wait();

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin();
          it != receiveProcessorComms.end(); ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      receiveProcessorComms.clear();

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end();
          ++it)
      {
        MPI_Type_free(&it->second.Type);
      }
      sendProcessorComms.clear();
      sendReceivePrepped = false;

      MPI_Waitall((int) (sendProcessorComms.size() + receiveProcessorComms.size()), &mRequests[0], &mStatuses[0]);

    }
  }
}
