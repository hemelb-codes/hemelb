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
          MPI_Gather(req->Pointer, 1, req->Type, NULL, 1, req->Type, it->first, communicator.GetCommunicator());
        }

      }
    }
  }
}

