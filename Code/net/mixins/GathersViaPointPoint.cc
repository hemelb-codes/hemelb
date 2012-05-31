#include "net/mixins/GathersViaPointPoint.h"
namespace hemelb
{
  namespace net
  {
    void GathersViaPointPoint::ReceiveGathers()
    {
      for (GatherProcComms::iterator receive_it = gatherReceiveProcessorComms.begin();
          receive_it != gatherReceiveProcessorComms.end(); ++receive_it)
      {
        int size;
        MPI_Type_size(receive_it->Type, &size);
        for (int source_rank = 0; source_rank <= communicator.GetSize(); source_rank++)
        {
          RequestReceive(&receive_it->Pointer + size * source_rank, 1, source_rank, receive_it->Type);
        }
      }
      gatherSendProcessorComms.clear();
    }
    void GathersViaPointPoint::SendGathers()
    {
      for (std::map<proc_t, GatherProcComms>::iterator send_it = gatherSendProcessorComms.begin();
          send_it != gatherSendProcessorComms.end(); ++send_it)
      {

        for (GatherProcComms::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
        {
          RequestSend(req->Pointer, 1, send_it->first, req->Type);
        }
      }

      gatherReceiveProcessorComms.clear();
    }
    void GathersViaPointPoint::ReceiveGatherVs()
    {
      for (GatherVReceiveProcComms::iterator receive_it = gatherVReceiveProcessorComms.begin();
          receive_it != gatherVReceiveProcessorComms.end(); ++receive_it)
      {
        for (int source_rank = 0; source_rank <= communicator.GetSize(); source_rank++)
        {
          RequestReceive(receive_it->Pointer + receive_it->Displacements[source_rank],
                         receive_it->Counts[source_rank],
                         source_rank,
                         receive_it->Type);
        }
      }

      gatherVReceiveProcessorComms.clear();
    }
    void GathersViaPointPoint::SendGatherVs()
    {
      for (std::map<proc_t, ProcComms>::iterator send_it = gatherVSendProcessorComms.begin();
          send_it != gatherVSendProcessorComms.end(); ++send_it)
      {

        for (ProcComms::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
        {
          RequestSend(req->Pointer, req->Count, send_it->first, req->Type);
        }

      }

      gatherVSendProcessorComms.clear();
    }
  }
}

