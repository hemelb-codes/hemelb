
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      HasCommsTestFixture::HasCommsTestFixture() {
	asyncCommQ = comm::Async::New(hemelbCommunicator);
      }

      void HasCommsTestFixture::Init(comm::Communicator::ConstPtr inst)
      {
	hemelbCommunicator = inst;
      }

      comm::Communicator::ConstPtr HasCommsTestFixture::Comms()
      {
	return hemelbCommunicator;
      }
      comm::Async::Ptr HasCommsTestFixture::Async()
      {
	return asyncCommQ;
      }
      comm::Communicator::ConstPtr HasCommsTestFixture::hemelbCommunicator = nullptr;
    }
  }
}

