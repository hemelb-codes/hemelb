
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/MockCommunicator.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {

      MockRequest::MockRequest() : comm::Request(), testCount(0) {
      }
      void MockRequest::Wait() {
      }
      bool MockRequest::Test() {
	if (testCount > 1) {
	  Wait();
	  return true;
	} else {
	  ++testCount;
	  return false;
	}
      }

      MockIbarrierRequest::MockIbarrierRequest(NullaryPredicate fn) : done(fn) {
      }
      void MockIbarrierRequest::Wait() {
	INFO("Predicate must be done before waiting!");
	REQUIRE(done());
      }
      bool MockIbarrierRequest::Test() {
	return done();
      }

      size_t MockRequestList::size() const {
	return reqs.size();
      }
      void MockRequestList::resize(size_t i) {
	reqs.resize(i);
      }
      void MockRequestList::push_back(Request::Ptr rp) {
	reqs.push_back(rp);
      }
	
      void MockRequestList::clear() {
	reqs.clear();
      }
	
      void MockRequestList::set(size_t i, Request::Ptr rp) {
	reqs[i] = rp;
      }
	
      void MockRequestList::WaitAll() {
	for (auto r: reqs) {
	  r->Wait();
	}
      }

      bool MockRequestList::TestAll() {
	bool alldone = true;
	for (auto r: reqs) {
	  if (!r->Test()) {
	    alldone = false;
	    break;
	  }
	}
	return alldone;
      }
            
      
      MockCommunicator::MockCommunicator(int rank_, int size_) :
	rank(rank_), size(size_)
      {

      }

      int MockCommunicator::Rank() const
      {
	return rank;
      }
      int MockCommunicator::Size() const
      {
	return size;
      }
	
#define NOTIMPLEMENTED FAIL("Communicator function not implemented")
	
      void MockCommunicator::Abort(int errCode) const { NOTIMPLEMENTED; }
	
      auto MockCommunicator::Duplicate() const -> Ptr { NOTIMPLEMENTED; return nullptr; }

      std::shared_ptr<Group> MockCommunicator::GetGroup() const { NOTIMPLEMENTED; return nullptr; }
      auto MockCommunicator::Create(std::shared_ptr<const Group> grp) const -> Ptr{ NOTIMPLEMENTED; return nullptr; }
      MpiFile MockCommunicator::OpenFile(const std::string& filename, int mode,
					 const MPI_Info info) const { NOTIMPLEMENTED; return MpiFile{}; }
      
      std::shared_ptr<RequestList> MockCommunicator::MakeRequestList() const {
	return std::make_shared<MockRequestList>();
      }
      void MockCommunicator::Barrier() const { NOTIMPLEMENTED; }
	
      std::shared_ptr<Request> MockCommunicator::Ibarrier() const {
	REQUIRE(ibarrier_conditions.size());
	auto& cond = ibarrier_conditions.front();
	auto ans = std::make_shared<MockIbarrierRequest>(cond);
	ibarrier_conditions.pop_front();
	return ans;
      }

      bool MockCommunicator::Iprobe(int source, int tag, MPI_Status* stat) const {
	// Look through the pending recvs and point to the first one that matches
	auto maybe_match = recv_data.begin();
	auto end = recv_data.end();
	for (; maybe_match != end; ++maybe_match) {
	  if (source == MPI_ANY_SOURCE || source == maybe_match->src) {
	    if (tag == MPI_ANY_TAG || tag == maybe_match->tag) {
	      break;
	    }
	  }
	}

	// No match - done
	if (maybe_match == end) {
	  return false;
	}
	  
	// Got a match, setup the status then move it to the front
	// of the queue
	stat->MPI_SOURCE = maybe_match->src;
	stat->MPI_TAG = maybe_match->tag;
	auto tmp = std::move(*maybe_match);
	recv_data.erase(maybe_match);
	recv_data.push_front(std::move(tmp));
	return true;
      }

      void MockCommunicator::ExpectationsAllCompleted() {
	REQUIRE(size_t(0) == gather_data.size());
	REQUIRE(size_t(0) == gatherv_data.size());
	REQUIRE(size_t(0) == send_data.size());
	REQUIRE(size_t(0) == ssend_data.size());
	REQUIRE(size_t(0) == recv_data.size());
	REQUIRE(size_t(0) == ibarrier_conditions.size());
      }
	

      void MockCommunicator::RequireIbarrier(MockIbarrierRequest::NullaryPredicate cond) {
	ibarrier_conditions.push_back(cond);
      }
      
      void MockCommunicator::BcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; }
      std::shared_ptr<Request> MockCommunicator::IbcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; return nullptr; }
      void MockCommunicator::AllreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; }
      std::shared_ptr<Request> MockCommunicator::IallreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; return nullptr; }
      std::shared_ptr<Request> MockCommunicator::IreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; return nullptr; }
      void MockCommunicator::ReduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; }
	
      void MockCommunicator::GatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
					void* recv, int recvcount, MPI_Datatype recvtype,
					int root) const {
	REQUIRE(gather_data.size());
	  
	auto& gd = gather_data.front();
	REQUIRE(gd.sendcount == sendcount);
	REQUIRE(gd.sendtype == sendtype);
	REQUIRE(gd.root == root);
	int elSize;
	MPI_Type_size(gd.sendtype, &elSize);
	  
	{ // assert that send[:] == relevant bit of stored data
	  int size_bytes = sendcount * elSize;
	  const char* expected = static_cast<const char*>(gd.send_buf);
	  const char* actual = static_cast<const char*>(send);
	  for (auto i = 0; i < size_bytes; ++i)
	    REQUIRE(expected[i] == actual[i]);
	}
	  
	if (Rank() == root) {
	  REQUIRE(gd.recvcount == recvcount);
	  REQUIRE(gd.recvtype == recvtype);

	  std::memcpy(recv, gd.ans_buf, recvcount * elSize * Size());
	}
	gather_data.pop_front();
      }
	
      void MockCommunicator::GathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
					 void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype,
					 int root) const {
	REQUIRE(gatherv_data.size());

	auto& gd = gatherv_data.front();
	REQUIRE(gd.sendcount == sendcount);
	REQUIRE(gd.sendtype == sendtype);
	REQUIRE(gd.root == root);
	  
	int elSize;
	MPI_Type_size(gd.sendtype, &elSize);
	  
	{ // assert that send[:] == relevant bit of stored data
	  int size_bytes = sendcount * elSize;
	  const char* expected = static_cast<const char*>(gd.send_buf);
	  const char* actual = static_cast<const char*>(sendbuf);
	  for (auto i = 0; i < size_bytes; ++i)
	    REQUIRE(expected[i] == actual[i]);
	}
	  
	if (Rank() == root) {
	  // The recv counts and displacements must match on root
	  for(auto i = 0; i < Size(); ++i) {
	    REQUIRE(gd.recvcounts[i] == recvcounts[i]);
	    REQUIRE(gd.displs[i] == displs[i]);
	  }
	  REQUIRE(gd.recvtype == recvtype);
	    
	  int nElTotal = std::accumulate(recvcounts, recvcounts+Size(), 0);
	  std::memcpy(recvbuf, gd.ans_buf, nElTotal*elSize);
	}
	gatherv_data.pop_front();
      }
	
      void MockCommunicator::AllgatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
					   void* recv, int recvcount, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
      void MockCommunicator::AllgathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
					    void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
      void MockCommunicator::ScatterImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
					 void* recvbuf, int recvcount, MPI_Datatype recvtype, int root) const { NOTIMPLEMENTED; }
      void MockCommunicator::AlltoallImpl(const void* send, int sendcount, MPI_Datatype sendtype,
					  void* recv, int recvcount, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
      void MockCommunicator::SendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				      int dest, int tag) const {
	REQUIRE(send_data.size());
	  
	auto& sd = send_data.front();
	REQUIRE(sd.count == sendcount);
	REQUIRE(sd.type == sendtype);
	REQUIRE(sd.dest == dest);
	REQUIRE(sd.tag == tag);
	send_data.pop_front();
      }
      void MockCommunicator::SsendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				       int dest, int tag) const {
	REQUIRE(ssend_data.size());
	  
	auto& sd = ssend_data.front();
	REQUIRE(sd.count == sendcount);
	REQUIRE(sd.type == sendtype);
	REQUIRE(sd.dest == dest);
	REQUIRE(sd.tag == tag);
	ssend_data.pop_front();
      }
	
      void MockCommunicator::RecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
				      int src, int tag, MPI_Status* stat) const {
	REQUIRE(recv_data.size());
	  
	auto& rd = recv_data.front();
	REQUIRE(rd.count == recvcount);
	REQUIRE(rd.type == recvtype);
	REQUIRE(rd.src == src);
	REQUIRE(rd.tag == tag);
	  
	int elSize;
	MPI_Type_size(rd.type, &elSize);
	  
	std::memcpy(recvbuf, rd.ans_buf, recvcount * elSize);
	recv_data.pop_front();

      }
	
      std::shared_ptr<Request> MockCommunicator::IsendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
							   int dest, int tag) const {
	SendImpl(sendbuf, sendcount, sendtype, dest, tag);
	return std::make_shared<MockRequest>();
      }
	
      std::shared_ptr<Request> MockCommunicator::IssendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
							    int dest, int tag) const {
	SsendImpl(sendbuf, sendcount, sendtype, dest, tag);
	return std::make_shared<MockRequest>();
      }
      std::shared_ptr<Request> MockCommunicator::IrecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
							   int source, int tag) const {
	RecvImpl(recvbuf, recvcount, recvtype, source, tag, MPI_STATUS_IGNORE);
	return std::make_shared<MockRequest>();
      }
      

    }
  }
}
