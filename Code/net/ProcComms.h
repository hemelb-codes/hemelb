
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_PROCCOMMS_H
#define HEMELB_NET_PROCCOMMS_H
#include "constants.h"
#include "net/mpi.h"
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

    class ProcComms : public BaseProcComms<SimpleRequest>
    {
      public:
        void CreateMPIType();
    };

    class GatherProcComms : public BaseProcComms<ScalarRequest>
    {

    };

    class AllToAllProcComms : public BaseProcComms<SimpleRequest>
    // Rank of request is not used - is all-to-all
    {

    };

    class GatherVReceiveProcComms : public BaseProcComms<GatherVReceiveRequest>
    {

    };
  }
}
#endif
