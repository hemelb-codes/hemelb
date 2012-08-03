// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_PHASED_STEPS_H
#define HEMELB_NET_PHASED_STEPS_H

namespace hemelb
{
  namespace net
  {
    namespace phased
    {
      namespace steps
      {
        enum Step
        {
          // Order significant here
          // BeginPhase must begin and EndPhase must end, those steps which should be called for a given phase.
          BeginAll = -1, // Called only before first phase
          BeginPhase = 0,
          Receive = 1,
          PreSend = 2,
          Send = 3,
          PreWait = 4,
          Wait = 5,
          EndPhase = 6,
          EndAll = 7, // Called only after final phase
        };
      }
    }
  }
}
#endif
