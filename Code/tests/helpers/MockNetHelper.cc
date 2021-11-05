// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/MockNetHelper.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {

      MockMpiCommunicator::MockMpiCommunicator(int rank_, int size_) :
	MpiCommunicator(rank_, size_)
      {
      }

      void MockNetHelper::setUp(const proc_t core_count, const proc_t current_core) {
	communicatorMock = new MockMpiCommunicator(current_core, core_count);
	netMock = new net::NetMock(*communicatorMock);
      }
      MockNetHelper::~MockNetHelper(){
	delete communicatorMock;
	delete netMock;
      }

    }
  }
}
