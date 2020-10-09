
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/MockCommsHelper.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      void MockCommsHelper::setUp(const proc_t core_count, const proc_t current_core)
      {
	communicator = std::make_shared<MockCommunicator>(current_core, core_count);
	commQ = comm::Async::New(communicator);
      }
	
      std::shared_ptr<MockCommunicator> MockCommsHelper::MockComms() {
	auto ans = std::dynamic_pointer_cast<MockCommunicator>(communicator);
	INFO("Cannot cast Communicator to a MockCommunicator");
	REQUIRE(ans);
	return ans;
      }
    }
  }
}

