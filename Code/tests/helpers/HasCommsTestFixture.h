
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H

#include "comm/Async.h"
#include "comm/Communicator.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      class HasCommsTestFixture
      {
      private:
	static comm::Communicator::ConstPtr hemelbCommunicator;
	comm::Async::Ptr asyncCommQ;
      public:
	static void Init(comm::Communicator::ConstPtr inst);
	HasCommsTestFixture();
      protected:
	static comm::Communicator::ConstPtr Comms();
	comm::Async::Ptr Async();
      };
    }
  }
}
#endif // HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H
