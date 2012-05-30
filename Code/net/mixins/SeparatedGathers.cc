#include "net/mixins/SeparatedGathers.h"
namespace hemelb
{
  namespace net
  {
    void SeparatedGathers::WaitGathers()
    {
      int gather_index = 0;
      for (std::vector<ScalarRequest>::iterator it = gatherReceiveProcessorComms.begin();
          it != gatherReceiveProcessorComms.end(); ++it)
      {
        ScalarRequest toself = gatherSendProcessorComms[communicator.GetRank()][gather_index];
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending/receiving Gather at core %i", communicator.GetRank());

        MPI_Gather(toself.Pointer,
                   1,
                   toself.Type,
                   it->Pointer,
                   1,
                   it->Type,
                   communicator.GetRank(),
                   communicator.GetCommunicator());
        ++gather_index;
      }

      for (std::map<proc_t, GatherProcComms>::iterator it = gatherSendProcessorComms.begin();
          it != gatherSendProcessorComms.end(); ++it)
      {
        if (it->first == communicator.GetRank())
        {
          continue;
        }
        for (std::vector<ScalarRequest>::iterator req = it->second.begin(); req != it->second.end(); req++)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Sending Gather at core %i to core %i",
                                                        communicator.GetRank(),
                                                        it->first);
          MPI_Gather(req->Pointer, 1, req->Type, NULL, 1, req->Type, it->first, communicator.GetCommunicator());
        }

      }
      gatherSendProcessorComms.clear();
      gatherReceiveProcessorComms.clear();
    }

    void SeparatedGathers::WaitGatherVs()
    {
      int gather_index = 0;
      for (std::vector<GatherVReceiveRequest>::iterator it = gatherVReceiveProcessorComms.begin();
          it != gatherVReceiveProcessorComms.end(); ++it)
      {
        BaseRequest toself = gatherVSendProcessorComms[communicator.GetRank()][gather_index];
        log::Logger::Log<log::Debug, log::OnePerCore>("Sending/receiving GatherV at core %i, getting %i, sending %i",
                                                      communicator.GetRank(),
                                                      it->Counts[0],
                                                      toself.Count);
        MPI_Gatherv(toself.Pointer,
                    toself.Count,
                    toself.Type,
                    it->Pointer,
                    it->Counts,
                    it->Displacements,
                    it->Type,
                    communicator.GetRank(),
                    communicator.GetCommunicator());
        ++gather_index;
      }

      for (std::map<proc_t, ProcComms>::iterator it = gatherVSendProcessorComms.begin();
          it != gatherVSendProcessorComms.end(); ++it)
      {
        if (it->first == communicator.GetRank())
        {
          continue;
        }
        for (std::vector<BaseRequest>::iterator req = it->second.begin(); req != it->second.end(); req++)
        {
          MPI_Gatherv(req->Pointer,
                      req->Count,
                      req->Type,
                      NULL,
                      NULL,
                      NULL,
                      req->Type,
                      it->first,
                      communicator.GetCommunicator());
        }

      }
      gatherVSendProcessorComms.clear();
      gatherVReceiveProcessorComms.clear();
    }
  }
}

