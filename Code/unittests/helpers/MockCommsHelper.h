
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKCOMMSHELPER_H
#define HEMELB_UNITTESTS_HELPERS_MOCKCOMMSHELPER_H

#include "unittests/helpers/MockCommunicator.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class MockCommsHelper
      {
      protected:
	MockCommsHelper()
	{
	}
	void setUp(const proc_t core_count, const proc_t current_core)
	{
	  communicator = std::make_shared<MockCommunicator>(current_core, core_count);
	  commQ = comm::Async::New(communicator);
	  //communicator = communicatorMock;
	}
	void tearDown()
	{
	}
	
	std::shared_ptr<MockCommunicator> MockComms() {
	  auto ans = std::dynamic_pointer_cast<MockCommunicator>(communicator);
	  CPPUNIT_ASSERT_MESSAGE("Cannot cast Communicator to a MockCommunicator", ans);
	  return ans;
	}
	//std::shared_ptr<MockCommunicator> communicatorMock;
	comm::Communicator::Ptr communicator;
	comm::Async::Ptr commQ;
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_MOCKCOMMSHELPER_H
