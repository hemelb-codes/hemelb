//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "timestep/CollectiveActor.h"
#include "Exception.h"

namespace hemelb
{
  namespace timestep
  {
    CollectiveActor::CollectiveActor(comm::Communicator::ConstPtr comm, reporting::Timer& wTimer) :
      isCollectiveRunning(false), mustWait(false), collectiveComm(comm->Duplicate()), waitTimer(wTimer), collectiveReq()
    {
    }
    
    void CollectiveActor::Receive()
    {
    }
    
    /**
     * Wait on the collectives to finish.
     */
    void CollectiveActor::Wait(void)
    {
      if (!isCollectiveRunning)
	throw Exception() << "Can only wait if the collective is running";
      
      waitTimer.Start();
      if (mustWait) {
        collectiveReq->Wait();
        mustWait = false;
        isCollectiveRunning = false;
      } else {
        bool done = collectiveReq->Test();
        if (done)
          isCollectiveRunning = false;
      }
      waitTimer.Stop();
    }
  }
}
