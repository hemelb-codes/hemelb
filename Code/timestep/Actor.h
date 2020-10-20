// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TIMESTEP_ACTOR_H
#define HEMELB_TIMESTEP_ACTOR_H

namespace hemelb
{
  namespace timestep
  {

    // An actor is a thing that performs its actions during a phase of
    // a timestep. See TimeStepManager.h for more details.
    class Actor
    {
    public:
      virtual ~Actor() = default;

    private:
      friend class TimeStepManager;

      virtual void BeginAll() = 0;
      virtual void Begin() = 0;
      virtual void Receive() = 0;
      virtual void PreSend() = 0;
      virtual void Send() = 0;
      virtual void PreWait() = 0;
      virtual void Wait() = 0;
      virtual void End() = 0;
      virtual void EndAll() = 0;
      
    };
  }
}

#endif
