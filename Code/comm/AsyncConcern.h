// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_ASYNC_CONCERN_H
#define HEMELB_COMM_ASYNC_CONCERN_H

#include "comm/Async.h"
#include "net/phased/Concern.h"
#include "net/phased/steps.h"

namespace hemelb
{
  namespace comm
  {
    using namespace net::phased;
    class AsyncConcern : public Concern
    {
    public:
      AsyncConcern(Async::Ptr async) : mAsync(async)
      {
      }
      bool CallAction(int action)
      {
	switch (static_cast<steps::Step>(action))
	{
	case steps::Wait:
	  mAsync->Wait();
	  return true;
	  
	default:
	  return false;
	}
      }
    private:
      Async::Ptr mAsync;
    };
    
  }
}

#endif
