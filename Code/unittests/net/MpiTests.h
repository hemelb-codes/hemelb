
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_NET_MPITESTS_H
#define HEMELB_UNITTESTS_NET_MPITESTS_H

#include <cppunit/TestFixture.h>
#include "net/mpi.h"

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      using namespace hemelb::net;
      /**
       * Unittests of the MPI abstraction layer
       *
       * This is hard to test with a single task....
       */
      class MpiTests : public CppUnit::TestFixture
      {
        public:
        CPPUNIT_TEST_SUITE (MpiTests);
        CPPUNIT_TEST (TestMpiComm);
        CPPUNIT_TEST_SUITE_END();

          void TestMpiComm()
          {
            MpiCommunicator commNull;
            CPPUNIT_ASSERT(!commNull);

            MpiCommunicator commWorld = MpiCommunicator::World();
            CPPUNIT_ASSERT(commWorld);

            CPPUNIT_ASSERT(commNull== commNull);
            CPPUNIT_ASSERT(commWorld == commWorld);
            CPPUNIT_ASSERT(commWorld != commNull);

            {
              MpiCommunicator commWorld2 = commWorld;
              CPPUNIT_ASSERT(commWorld2 == commWorld);
            }

            {
              MpiCommunicator commWorld2 = MpiCommunicator::World();
              CPPUNIT_ASSERT(commWorld2 == commWorld);
            }
            {
              MpiGroup groupWorld = commWorld.Group();
              MpiCommunicator commWorld2 = commWorld.Create(groupWorld);
              // Same ranks, but different context.
              CPPUNIT_ASSERT(commWorld2 != commWorld);
            }
          }
      };
      CPPUNIT_TEST_SUITE_REGISTRATION (MpiTests);
    }
  }
}
#endif // HEMELB_UNITTESTS_NET_MPITESTS_H
