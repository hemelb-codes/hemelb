
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "tests/helpers/MockCommunicator.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      TEST_CASE("MockCommunicatorTests") {

	SECTION("Simple") {
	  auto comm = std::make_shared<MockCommunicator>(0, 1023);
	  REQUIRE(0 == comm->Rank());
	  REQUIRE(1023 == comm->Size());
	}

	SECTION("Send") {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(0, 2);
	  int number = 42;
	  comm->RequireSend(number, 1, 0);
	  comm->Send(42, 1);
	  comm->ExpectationsAllCompleted();
	}

	SECTION("Isend") {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(0, 2);
	  int number = 42;
	  comm->RequireSend(number, 1, 0);
	  auto req = comm->Isend(42, 1);
	  req->Wait();
	  comm->ExpectationsAllCompleted();
	}

	SECTION("Recv") {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(1, 2);
	  int number = 42;
	  comm->RequireRecv(number, 1, 0);
	  int ans;
	  comm->Recv(ans, 1);
	  comm->ExpectationsAllCompleted();
	  REQUIRE(number == ans);
	}

	SECTION("Irecv") {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(1, 2);
	  int number = 42;
	  comm->RequireRecv(number, 1, 0);
	  int ans;
	  auto req = comm->Irecv(ans, 1);
	  req->Wait();
	  comm->ExpectationsAllCompleted();
	  REQUIRE(number == ans);
	}

	SECTION("GatherRoot") {
	  const int rank = 2;
	  const int size = 8;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  std::vector<int> vals{9,8,7,6,5,4,3,2};
	  comm->RequireGather(vals, 2);

	  int myval = 7;
	  auto ans = comm->Gather(myval, 2);
	  
	  REQUIRE(size == int(ans.size()));
	  comm->ExpectationsAllCompleted();
	}
	SECTION("GatherNonRoot") {
	  const int rank = 0;
	  const int size = 8;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  std::vector<int> vals{9,8,7,6,5,4,3,2};
	  comm->RequireGather(vals, 2);

	  int myval = 9;
	  auto ans = comm->Gather(myval, 2);
	  
	  REQUIRE(0 == int(ans.size()));
	  comm->ExpectationsAllCompleted();
	}

	SECTION("GatherV") {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);

	  std::vector<int> counts{1,2,3,0};
	  std::vector<char> data{3, 4,5, 6,7,8};
	  comm->RequireGatherV(data, counts, 0);

	  std::vector<char> mydata{3};
	  auto ans = comm->GatherV(mydata, counts, 0);
	  REQUIRE(data == ans);
	  comm->ExpectationsAllCompleted();
	}

	SECTION("Iprobe") {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);

	  int val1 = 42;
	  int val2 = 24;
	  comm->RequireRecv(val1, 1, 0);
	  comm->RequireRecv(val2, 2, 0);
	  // message q has recvs ordered (1,2)
	  MPI_Status stat;
	  REQUIRE(comm->Iprobe(2, MPI_ANY_TAG, &stat));
	  // messages should now be ordered (2,1)
	  int ans;
	  comm->Recv(ans, 2);
	  REQUIRE(val2 == ans);
	  
	  comm->Recv(ans, 1);
	  REQUIRE(val1 == ans);
	  
	  REQUIRE(!comm->Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, &stat));
	  comm->ExpectationsAllCompleted();
	}

	SECTION("Ibarrier") {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  // Try a dead simple one first
	  comm->RequireIbarrier([]() {
	      return true;
	    });
	  auto req = comm->Ibarrier();
	  REQUIRE(req->Test());

	  // Now this guy will count the calls and only work on the 3rd go
	  int i = 0;
	  comm->RequireIbarrier([&i]() {
	      return ++i == 3;
	    });
	  req = comm->Ibarrier();
	  REQUIRE(!req->Test());
	  REQUIRE(!req->Test());
	  req->Wait();
	  
	  comm->ExpectationsAllCompleted();
	}
      }
    }
  }
}

