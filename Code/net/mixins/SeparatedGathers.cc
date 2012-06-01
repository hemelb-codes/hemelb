#include "net/mixins/SeparatedGathers.h"
namespace hemelb
{
  namespace net
  {
    void SeparatedGathers::WaitGathers()
    {

      for (std::map<proc_t, GatherProcComms>::iterator send_it = gatherSendProcessorComms.begin();
          send_it != gatherSendProcessorComms.end(); ++send_it)
      {
        if (send_it->first == communicator.GetRank())
        {
          int gather_index = 0;
          for (GatherProcComms::iterator receive_it = gatherReceiveProcessorComms.begin();
              receive_it != gatherReceiveProcessorComms.end(); ++receive_it)
          {
            ScalarRequest toself = send_it->second[gather_index];

            MPI_Gather(toself.Pointer,
                       1,
                       toself.Type,
                       receive_it->Pointer,
                       1,
                       receive_it->Type,
                       communicator.GetRank(),
                       communicator.GetCommunicator());
            ++gather_index;
          }
        }
        else
        {
          for (GatherProcComms::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
          {
            MPI_Gather(req->Pointer, 1, req->Type, NULL, 1, req->Type, send_it->first, communicator.GetCommunicator());
          }
        }
      }
      gatherSendProcessorComms.clear();
      gatherReceiveProcessorComms.clear();
    }

    void SeparatedGathers::WaitGatherVs()
    {

      for (std::map<proc_t, ProcComms>::iterator send_it = gatherVSendProcessorComms.begin();
          send_it != gatherVSendProcessorComms.end(); ++send_it)
      {
        if (send_it->first == communicator.GetRank())
        {
          int gather_index = 0;
          for (GatherVReceiveProcComms::iterator receive_it = gatherVReceiveProcessorComms.begin();
              receive_it != gatherVReceiveProcessorComms.end(); ++receive_it)
          {
            BaseRequest toself = send_it->second[gather_index];

            MPI_Gatherv(toself.Pointer,
                        toself.Count,
                        toself.Type,
                        receive_it->Pointer,
                        receive_it->Counts,
                        receive_it->Displacements,
                        receive_it->Type,
                        communicator.GetRank(),
                        communicator.GetCommunicator());
            ++gather_index;
          }
        }
        else
        {

          for (ProcComms::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
          {
            MPI_Gatherv(req->Pointer,
                        req->Count,
                        req->Type,
                        NULL,
                        NULL,
                        NULL,
                        req->Type,
                        send_it->first,
                        communicator.GetCommunicator());
          }

        }
      }
      gatherVSendProcessorComms.clear();
      gatherVReceiveProcessorComms.clear();
    }
  }
}

