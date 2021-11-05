// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_NET_PHASED_MOCKITERATEDACTION_H
#define HEMELB_TESTS_NET_PHASED_MOCKITERATEDACTION_H

#include <string>
#include <sstream>

#include "net/IteratedAction.h"

namespace hemelb
{
  namespace tests
  {

    class MockIteratedAction : public hemelb::net::IteratedAction
    {
    public:
      MockIteratedAction(const std::string & name);

      void PreSend();
      void PreReceive();
      void PostReceive();
      void EndIteration();
      void RequestComms();
      void Reset();
      std::string CallsSoFar();

    private:
      std::string name;
      std::stringstream calls;
    };
  }
}

#endif
