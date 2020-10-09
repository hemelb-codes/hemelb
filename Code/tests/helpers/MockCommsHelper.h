
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_MOCKCOMMSHELPER_H
#define HEMELB_TESTS_HELPERS_MOCKCOMMSHELPER_H

#include "units.h"
#include "comm/Communicator.h"
#include "comm/Async.h"
#include "tests/helpers/MockCommunicator.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      class MockCommsHelper
      {
      protected:
	comm::Communicator::Ptr communicator = nullptr;
	comm::Async::Ptr commQ = nullptr;
	void setUp(const proc_t core_count, const proc_t current_core);
	std::shared_ptr<MockCommunicator> MockComms();
      };

    }
  }
}

#endif // HEMELB_TESTS_HELPERS_MOCKCOMMSHELPER_H
