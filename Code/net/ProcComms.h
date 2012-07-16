// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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
