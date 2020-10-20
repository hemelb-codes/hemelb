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
      inline AsyncConcern(Async::Ptr async) : mAsync(async)
      {
      }
      ~AsyncConcern() = default;

    private:
      inline void BeginAll() override {
      }
      inline void Begin() override {
      }
      inline void Receive() override {
      }
      inline void PreSend() override {
      }
      inline void Send() override {
      }
      inline void PreWait() override {
      }

      inline void Wait() override {
	mAsync->Wait();
      }

      inline void End() override {
      }
      inline void EndAll() override {
      }

      Async::Ptr mAsync;
    };
    
  }
}

#endif
