//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "net/CollectiveAction.h"
#include "net/phased/steps.h"

namespace hemelb
{
  namespace net
  {
    CollectiveAction::CollectiveAction(const MpiCommunicator& comm, reporting::Timer& wTimer) :
        collectiveComm(comm.Duplicate()), waitTimer(wTimer), collectiveReq()
    {
    }
    bool CollectiveAction::CallAction(int action)
    {
      switch (static_cast<phased::steps::Step>(action))
      {
        case phased::steps::BeginPhase:
          RequestComms();
          return true;
        case phased::steps::PreSend:
          PreSend();
          return true;
        case phased::steps::Send:
          Send();
          return true;
        case phased::steps::PreWait:
          PreReceive();
          return true;
        case phased::steps::Wait:
          Wait();
          return true;
        case phased::steps::EndPhase:
          PostReceive();
          return true;
        case phased::steps::EndAll:
          EndIteration();
          return true;
        default:
          return false;
      }
    }
    /**
     * Wait on the collectives to finish.
     */
    void CollectiveAction::Wait(void)
    {
      waitTimer.Start();
      collectiveReq.Wait();
      waitTimer.Stop();
    }

  }
}
