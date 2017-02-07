// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_ASYNC_CONCERN_H
#define HEMELB_COMM_ASYNC_CONCERN_H

#include "comm/Async.h"
#include "timestep/Actor.h"

namespace hemelb
{
  namespace comm
  {
    class AsyncConcern : public timestep::Actor
    {
    public:
      AsyncConcern(Async::Ptr async) : mAsync(async)
      {
      }
      inline virtual void BeginAll() {
      }
      inline virtual void Begin()  {
      }
      inline virtual void Receive()  {
      }
      inline virtual void PreSend()  {
      }
      inline virtual void Send()  {
      }
      inline virtual void PreWait()  {
      }
      
      inline virtual void Wait() {
	mAsync->Wait();
      }
      
      inline virtual void End()  {
      }
      inline virtual void EndAll()  {
      }
      
    private:
      Async::Ptr mAsync;
    };
    
  }
}

#endif
