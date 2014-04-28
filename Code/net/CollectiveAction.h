//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_COLLECTIVEACTION_H
#define HEMELB_NET_COLLECTIVEACTION_H

#include "net/IteratedAction.h"
#include "reporting/Timers.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace net
  {
    class CollectiveAction : public IteratedAction
    {
      public:
        bool CallAction(int action);

        /**
         * Initiate the collective.
         */
        virtual void Send(void) = 0;

        /**
         * Wait on the collectives to finish.
         */
        virtual void Wait(void);

      protected:
        CollectiveAction(const MpiCommunicator& comm, reporting::Timers& timings);

        /**
         * Private communicator for non-blocking collectives.
         */
        MpiCommunicator collectiveComm;
        /**
         * Timings for the wait etc.
         */
        reporting::Timers& timings;
        /**
         * Request object for the collective
         */
        MpiRequest collectiveReq;
     };
  }
}

#endif /* HEMELB_NET_COLLECTIVEACTION_H */
