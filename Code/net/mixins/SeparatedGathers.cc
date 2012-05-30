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
          for (std::vector<ScalarRequest>::iterator receive_it = gatherReceiveProcessorComms.begin();
              receive_it != gatherReceiveProcessorComms.end(); ++receive_it)
          {
            ScalarRequest toself = send_it->second[gather_index];
            log::Logger::Log<log::Debug, log::OnePerCore>("Sending/receiving Gather at core %i",
                                                          communicator.GetRank());

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
          for (std::vector<ScalarRequest>::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Sending Gather at core %i to core %i",
                                                          communicator.GetRank(),
                                                          send_it->first);
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
          for (std::vector<GatherVReceiveRequest>::iterator receive_it = gatherVReceiveProcessorComms.begin();
              receive_it != gatherVReceiveProcessorComms.end(); ++receive_it)
          {
            BaseRequest toself = send_it->second[gather_index];
            log::Logger::Log<log::Debug, log::OnePerCore>("Sending/receiving GatherV at core %i, getting %i, sending %i",
                                                          communicator.GetRank(),
                                                          receive_it->Counts[0],
                                                          toself.Count);
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

          for (std::vector<BaseRequest>::iterator req = send_it->second.begin(); req != send_it->second.end(); req++)
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

