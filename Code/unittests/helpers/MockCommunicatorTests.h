
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKCOMMUNICATORTESTS_H
#define HEMELB_UNITTESTS_HELPERS_MOCKCOMMUNICATORTESTS_H

#include <cppunit/TestFixture.h>
#include "unittests/helpers/MockCommunicator.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      class MockCommunicatorTests : public CppUnit::TestFixture
      {
	CPPUNIT_TEST_SUITE( MockCommunicatorTests);
	CPPUNIT_TEST( TestSimple);
	
	CPPUNIT_TEST( TestSend);
	CPPUNIT_TEST( TestIsend);
	
	CPPUNIT_TEST( TestRecv);
	CPPUNIT_TEST( TestIrecv);

	CPPUNIT_TEST( TestGatherRoot);
	CPPUNIT_TEST( TestGatherNonRoot);

	CPPUNIT_TEST( TestGatherV);

	CPPUNIT_TEST( TestIprobe);
	CPPUNIT_TEST( TestIbarrier);
	
	CPPUNIT_TEST_SUITE_END();
	
        public:
	
	void TestSimple() {
	  auto comm = std::make_shared<MockCommunicator>(0, 1023);
	  CPPUNIT_ASSERT_EQUAL(0, comm->Rank());
	  CPPUNIT_ASSERT_EQUAL(1023, comm->Size());
	}

	void TestSend() {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(0, 2);
	  int number = 42;
	  comm->RequireSend(number, 1, 0);
	  comm->Send(42, 1);
	  comm->ExpectationsAllCompleted();
	}

	void TestIsend() {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(0, 2);
	  int number = 42;
	  comm->RequireSend(number, 1, 0);
	  auto req = comm->Isend(42, 1);
	  req->Wait();
	  comm->ExpectationsAllCompleted();
	}

	void TestRecv() {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(1, 2);
	  int number = 42;
	  comm->RequireRecv(number, 1, 0);
	  int ans;
	  comm->Recv(ans, 1);
	  comm->ExpectationsAllCompleted();
	  CPPUNIT_ASSERT_EQUAL(number, ans);
	}

	void TestIrecv() {
	  // Send a message 0 -> 1
	  auto comm = std::make_shared<MockCommunicator>(1, 2);
	  int number = 42;
	  comm->RequireRecv(number, 1, 0);
	  int ans;
	  auto req = comm->Irecv(ans, 1);
	  req->Wait();
	  comm->ExpectationsAllCompleted();
	  CPPUNIT_ASSERT_EQUAL(number, ans);
	}

	void TestGatherRoot() {
	  const int rank = 2;
	  const int size = 8;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  std::vector<int> vals{9,8,7,6,5,4,3,2};
	  comm->RequireGather(vals, 2);

	  int myval = 7;
	  auto ans = comm->Gather(myval, 2);
	  
	  CPPUNIT_ASSERT_EQUAL(size, int(ans.size()));
	  comm->ExpectationsAllCompleted();
	}
	void TestGatherNonRoot() {
	  const int rank = 0;
	  const int size = 8;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  std::vector<int> vals{9,8,7,6,5,4,3,2};
	  comm->RequireGather(vals, 2);

	  int myval = 9;
	  auto ans = comm->Gather(myval, 2);
	  
	  CPPUNIT_ASSERT_EQUAL(0, int(ans.size()));
	  comm->ExpectationsAllCompleted();
	}

	void TestGatherV() {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);

	  std::vector<int> counts{1,2,3,0};
	  std::vector<char> data{3, 4,5, 6,7,8};
	  comm->RequireGatherV(data, counts, 0);

	  std::vector<char> mydata{3};
	  auto ans = comm->GatherV(mydata, counts, 0);
	  CPPUNIT_ASSERT_EQUAL(data, ans);
	  comm->ExpectationsAllCompleted();
	}

	void TestIprobe() {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);

	  int val1 = 42;
	  int val2 = 24;
	  comm->RequireRecv(val1, 1, 0);
	  comm->RequireRecv(val2, 2, 0);
	  // message q has recvs ordered (1,2)
	  MPI_Status stat;
	  CPPUNIT_ASSERT(comm->Iprobe(2, MPI_ANY_TAG, &stat));
	  // messages should now be ordered (2,1)
	  int ans;
	  comm->Recv(ans, 2);
	  CPPUNIT_ASSERT_EQUAL(val2, ans);
	  
	  comm->Recv(ans, 1);
	  CPPUNIT_ASSERT_EQUAL(val1, ans);
	  
	  CPPUNIT_ASSERT(!comm->Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, &stat));
	  comm->ExpectationsAllCompleted();
	}

	void TestIbarrier() {
	  const int rank = 0;
	  const int size = 4;
	  auto comm = std::make_shared<MockCommunicator>(rank, size);
	  
	  // Try a dead simple one first
	  comm->RequireIbarrier([]() {
	      return true;
	    });
	  auto req = comm->Ibarrier();
	  CPPUNIT_ASSERT(req->Test());

	  // Now this guy will count the calls and only work on the 3rd go
	  int i = 0;
	  comm->RequireIbarrier([&i]() {
	      return ++i == 3;
	    });
	  req = comm->Ibarrier();
	  CPPUNIT_ASSERT(!req->Test());
	  CPPUNIT_ASSERT(!req->Test());
	  req->Wait();
	  
	  comm->ExpectationsAllCompleted();
	}
      };
      CPPUNIT_TEST_SUITE_REGISTRATION( MockCommunicatorTests);
    }
  }
}

#endif
