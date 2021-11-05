// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_NET_PHASED_MOCKCONCERN_H
#define HEMELB_TESTS_NET_PHASED_MOCKCONCERN_H

#include <string>
#include <vector>

#include "net/phased/Concern.h"

namespace hemelb
{
  namespace tests
  {
    class MockConcern : public net::phased::Concern
    {
    public:
      MockConcern(const std::string &name);

      bool CallAction(int action);

      std::vector<int> const & ActionsCalled() const;

    private:
      std::vector<int> calls;
      std::string name;
    };

  }
}
#endif
