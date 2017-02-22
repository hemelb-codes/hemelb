// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TIMESTEP_COLLECTIVEACTOR_H
#define HEMELB_TIMESTEP_COLLECTIVEACTOR_H

#include "timestep/Actor.h"
#include "reporting/Timers.h"
#include "comm/Request.h"

namespace hemelb
{
  namespace timestep
  {
    class CollectiveActor : public Actor
    {
    public:
      inline void MustFinishThisTimeStep()
      {

	mustWait = true;
      }

    protected:
      CollectiveActor(comm::Communicator::ConstPtr comm, reporting::Timer& waitTimer);
      
      // This should be a no-op for a collective
      virtual void Receive();
      
      // Subclass must supply the send method with the collective
      // virtual void Send();
      
      // Progress the communication
      virtual void Wait(void) final;

      bool isCollectiveRunning;
      bool mustWait;

      /**
       * Private communicator for non-blocking collectives.
       */
      comm::Communicator::Ptr collectiveComm;
      /**
       * Timings for the wait etc.
       */
      reporting::Timer& waitTimer;
      /**
       * Request object for the collective
       */
      comm::Request::Ptr collectiveReq;
    };

  }
}

#endif /* HEMELB_NET_COLLECTIVEACTION_H */
