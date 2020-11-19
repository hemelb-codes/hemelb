// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/net/phased/MockIteratedAction.h"
namespace hemelb
{
  namespace tests
  {

    MockIteratedAction::MockIteratedAction(const std::string & name) :
      name(name), calls()
    {
    }
    void MockIteratedAction::PreSend()
    {
      calls << "PreSend, " << std::flush;
    }
    void MockIteratedAction::PreReceive()
    {
      calls << "PreReceive, " << std::flush;
    }
    void MockIteratedAction::PostReceive()
    {
      calls << "PostReceive, " << std::flush;
    }
    void MockIteratedAction::EndIteration()
    {
      calls << "EndIteration, " << std::flush;
    }
    void MockIteratedAction::RequestComms()
    {
      calls << "RequestComms, " << std::flush;
    }
    void MockIteratedAction::Reset()
    {
      calls << "Reset, " << std::flush;
    }
    std::string MockIteratedAction::CallsSoFar(){
      return calls.str();
    }

  }
}
