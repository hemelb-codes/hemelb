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
    CollectiveAction::CollectiveAction(comm::Communicator::ConstPtr comm, reporting::Timer& wTimer) :
        isCollectiveRunning(false), mustWait(false), collectiveComm(comm->Duplicate()), waitTimer(wTimer), collectiveReq()
    {
    }

    bool CollectiveAction::CallAction(int action)
    {
      if (isCollectiveRunning)
      {
        switch (static_cast<phased::steps::Step>(action))
        {
          case phased::steps::Wait:
            Wait();
            return true;
          default:
            return false;
        }
      }
      else
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
            isCollectiveRunning = true;
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
    }

    /**
     * Wait on the collectives to finish.
     */
    void CollectiveAction::Wait(void)
    {
      waitTimer.Start();
      if (mustWait) {
        collectiveReq->Wait();
        mustWait = false;
        isCollectiveRunning = false;
      }
      else
      {
        bool done = collectiveReq->Test();
        if (done)
          isCollectiveRunning = false;
        }
      waitTimer.Stop();
    }
  }
}
