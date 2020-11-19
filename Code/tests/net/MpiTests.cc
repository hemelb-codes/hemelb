// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "net/mpi.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::net;
    /**
     * Unittests of the MPI abstraction layer
     *
     * This is hard to test with a single task....
     */
    TEST_CASE("MpiTests") {
      MpiCommunicator commNull;
      REQUIRE(!commNull);

      MpiCommunicator commWorld = MpiCommunicator::World();
      REQUIRE(commWorld);

      REQUIRE(commNull== commNull);
      REQUIRE(commWorld == commWorld);
      REQUIRE(commWorld != commNull);

      SECTION("Copy assignment works") {
	MpiCommunicator commWorld2 = commWorld;
	REQUIRE(commWorld2 == commWorld);
      }

      SECTION("World factor returns objects that compare equal") {
	MpiCommunicator commWorld2 = MpiCommunicator::World();
	REQUIRE(commWorld2 == commWorld);
      }

      SECTION("Comms with the same group are different") {
	MpiGroup groupWorld = commWorld.Group();
	MpiCommunicator commWorld2 = commWorld.Create(groupWorld);
	// Same ranks, but different context.
	REQUIRE(commWorld2 != commWorld);
      }
    }
  }
}
