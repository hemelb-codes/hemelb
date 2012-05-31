#ifndef HEMELB_NET_PROCCOMMS_H
#define HEMELB_NET_PROCCOMMS_H
#include "constants.h"
#include "mpiInclude.h"
#include "net/StoredRequest.h"
#include <deque>
#include <vector>

namespace hemelb
{
  namespace net
  {
    template<class Request>
    class BaseProcComms : public std::deque<Request>
    {
      public:
        MPI_Datatype Type;
    };

    class ProcComms : public BaseProcComms<BaseRequest>
    {
      public:
        void CreateMPIType();
    };

    class GatherProcComms : public BaseProcComms<ScalarRequest>
    {

    };

    class GatherVReceiveProcComms : public BaseProcComms<GatherVReceiveRequest>
    {

    };
  }
}
#endif
