// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_MOCKNETHELPER_H
#define HEMELB_TESTS_HELPERS_MOCKNETHELPER_H

#include "net/MpiCommunicator.h"

#include "tests/net/NetMock.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      class MockMpiCommunicator : public net::MpiCommunicator
      {
      public:
	/***
	 * Constructor for a dummy communicator
	 * Can be useful for testing but can't actually be used
	 * @param rank
	 * @param size
	 */
	MockMpiCommunicator(int rank_, int size_);
      };

      class MockNetHelper
      {
      protected:
	net::MpiCommunicator* communicatorMock = nullptr;
	net::NetMock* netMock = nullptr;
	void setUp(const proc_t core_count, const proc_t current_core);
	~MockNetHelper();

      };

    }
  }
}

#endif // HEMELB_TESTS_HELPERS_RANDOMSOURCE_H
