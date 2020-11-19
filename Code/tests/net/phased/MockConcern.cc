// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/net/phased/MockConcern.h"

namespace hemelb
{
  namespace tests
  {

    MockConcern::MockConcern(const std::string &name) :
      calls(), name(name)
    {
    }

    bool MockConcern::CallAction(int action)
    {
      // this is where a real concern would switch on the action, and call the appropriate method
      calls.push_back(action);
      return true;
    }

    std::vector<int> const & MockConcern::ActionsCalled() const
    {
      return calls;
    }

  }
}
