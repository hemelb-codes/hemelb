// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "timestep/Actor.h"

namespace hemelb
{
  namespace timestep
  {
    void Actor::BeginAll() {}
    void Actor::Begin() {}
    void Actor::Receive() {}
    void Actor::PreSend() {}
    void Actor::Send() {}
    void Actor::PreWait() {}
    void Actor::Wait() {}
    void Actor::End() {}
    void Actor::EndAll() {}  
  }
}

