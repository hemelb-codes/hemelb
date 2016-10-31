
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKCOMMUNICATOR_H
#define HEMELB_UNITTESTS_HELPERS_MOCKCOMMUNICATOR_H

#include <numeric>
#include "comm/Communicator.h"
#include "comm/Request.h"

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      using namespace comm;

      class MockRequest : public comm::Request
      {
      public:
	MockRequest() : comm::Request(), testCount(0) {
	}
        virtual void Wait() {
	  
	}
        virtual bool Test() {
	  if (testCount > 1) {
	    Wait();
	    return true;
	  } else {
	    ++testCount;
	    return false;
	  }
	}
	
      private:
	int testCount;
      };

      class MockIbarrierRequest : public MockRequest
      {
      public:
	typedef std::function<bool()> NullaryPredicate;
	MockIbarrierRequest(NullaryPredicate fn) : done(fn){
	}
	virtual void Wait() {
	  CPPUNIT_ASSERT_MESSAGE("Predicate must be done before waiting!", done());
	}
	virtual bool Test() {
	  return done();
	}
      private:
	NullaryPredicate done;
      };
      
      class MockRequestList : public comm::RequestList
      {
	std::vector<Request::Ptr> reqs;
      public:
	virtual size_t size() const {
	  return reqs.size();
	}
	virtual void resize(size_t i) {
	  reqs.resize(i);
	}
	virtual void push_back(Request::Ptr rp) {
	  reqs.push_back(rp);
	}
	
	virtual void clear() {
	  reqs.clear();
	}
	
	virtual void set(size_t i, Request::Ptr rp) {
	  reqs[i] = rp;
	}
	
	virtual void WaitAll() {
	  for (auto r: reqs) {
	    r->Wait();
	  }
	}
	virtual bool TestAll() {
	  bool alldone = true;
	  for (auto r: reqs) {
	    if (!r->Test()) {
	      alldone = false;
	      break;
	    }
	  }
	  return alldone;
	}
      };
      
      // class MockSendReq : public MockRequest
      // {
      // 	MockSendReq(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
      // 		    int dest, int tag) : mSendbuf(sendbuf), mSendcount(sendcount), mSendtype(sendtype),
      // 					 mDest(dest), mTag(tag)
      // 	{
      // 	}
	
      // private:
      // 	const void *mSendbuf;
      // 	int mSendcount;
      // 	MPI_Datatype mSendtype;
      // 	int mDest;
      // 	int mTag;
      // };
      
      class MockCommunicator : public comm::Communicator
      {
      public:
	/***
	 * Constructor for a dummy communicator
	 * Can be useful for testing but can't actually be used
	 * @param rank
	 * @param size
	 */
	MockCommunicator(int rank_, int size_) :
	  rank(rank_), size(size_)
	{

	}

	virtual inline int Rank() const
	{
	  return rank;
	}
	virtual inline int Size() const
	{
	  return size;
	}
	
#define NOTIMPLEMENTED CPPUNIT_FAIL("Communicator function not implemented")
	
	virtual void Abort(int errCode) const { NOTIMPLEMENTED; }
	
        virtual Ptr Duplicate() const { NOTIMPLEMENTED; return nullptr; }

        virtual std::shared_ptr<Group> GetGroup() const { NOTIMPLEMENTED; return nullptr; }
        virtual Ptr Create(std::shared_ptr<const Group> grp) const { NOTIMPLEMENTED; return nullptr; }
        virtual std::shared_ptr<MpiFile> OpenFile(const std::string& filename, int mode,
						  const MPI_Info info = MPI_INFO_NULL) const { NOTIMPLEMENTED; return nullptr; }
      
        virtual std::shared_ptr<RequestList> MakeRequestList() const {
	  return std::make_shared<MockRequestList>();
	}
        virtual void Barrier() const { NOTIMPLEMENTED; }
	
        virtual std::shared_ptr<Request> Ibarrier() const {
	  CPPUNIT_ASSERT(ibarrier_conditions.size());
	  auto& cond = ibarrier_conditions.front();
	  auto ans = std::make_shared<MockIbarrierRequest>(cond);
	  ibarrier_conditions.pop_front();
	  return ans;
	}

        virtual bool Iprobe(int source, int tag, MPI_Status* stat=MPI_STATUS_IGNORE) const {
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

	void ExpectationsAllCompleted() {
	  CPPUNIT_ASSERT_EQUAL(size_t(0), gather_data.size());
	  CPPUNIT_ASSERT_EQUAL(size_t(0), gatherv_data.size());
	  CPPUNIT_ASSERT_EQUAL(size_t(0), send_data.size());
	  CPPUNIT_ASSERT_EQUAL(size_t(0), ssend_data.size());
	  CPPUNIT_ASSERT_EQUAL(size_t(0), recv_data.size());
	  CPPUNIT_ASSERT_EQUAL(size_t(0), ibarrier_conditions.size());
	}
	
	template<class T>
	void RequireGather(const std::vector<T>& data, int root)
	{
	  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of elements must be a multiple of comm size",
				       0, int(data.size()) % Size());
	  int sendcount = data.size() / Size();
	  gather_data.push_back({data.data() + sendcount*Rank(), sendcount, comm::MpiDataType<T>(),
		data.data(), sendcount, comm::MpiDataType<T>(),
		root});
	}
	
	template<class T>
	void RequireGatherV(const std::vector<T>& data, const std::vector<int>& recvcounts, int root)
	{
	  CPPUNIT_ASSERT_EQUAL(Size(), int(recvcounts.size()));
	  std::vector<int> displs(Size());
	  int total = 0;
	  for(size_t i = 0; i < Size(); ++i) {
	    displs[i] = total;
	    total += recvcounts[i];
	  }
	  auto my_send_data = data.data() + displs[Rank()];
	  CPPUNIT_ASSERT_EQUAL_MESSAGE("Size of data does not match sum of recvcounts",
				       size_t(total), data.size());
	  gatherv_data.push_back({my_send_data, recvcounts[Rank()], comm::MpiDataType<T>(),
		data.data(), recvcounts.data(), std::move(displs), comm::MpiDataType<T>(),
		root});
	}
	
	template <typename T>
	void RequireSend(const T* valPtr, int count, int dest, int tag)
	{
	  static_assert(!std::is_pointer<T>::value,
			"You are trying to require the sending of a pointer type");
	  static_assert(!is_std_vector<T>::value,
			"You are trying to require the sending of a std::vector type");

	  CPPUNIT_ASSERT(dest >= 0);
	  CPPUNIT_ASSERT(dest < Size());
	  send_data.push_back({valPtr, count, comm::MpiDataType<T>(), dest, tag});
	}
	template <typename T>
	void RequireSend(const T& val, int dest, int tag)
	{
	  RequireSend(&val, 1, dest, tag);
	}
	template <typename T>
	void RequireSend(const std::vector<T>& vals, int dest, int tag)
	{
	  RequireSend(vals.data(), vals.size(), dest, tag);
	}
	
	template <typename T>
	void RequireRecv(const T* valPtr, int count, int src, int tag)
	{
	  static_assert(!std::is_pointer<T>::value,
			"You are trying to require the receiving of a pointer type");
	  static_assert(!is_std_vector<T>::value,
			"You are trying to require the receiving of a std::vector type");
	  CPPUNIT_ASSERT(src >= 0);
	  CPPUNIT_ASSERT(src < Size());
	  recv_data.push_back({valPtr, count, comm::MpiDataType<T>(), src, tag});	  
	}
	template <typename T>
	void RequireRecv(const T& val, int src, int tag)
	{
	  RequireRecv(&val, 1, src, tag);
	}
	template <typename T>
	void RequireRecv(const std::vector<T>& vals, int src, int tag)
	{
	  RequireRecv(vals.data(), vals.size(), src, tag);
	}

	void RequireIbarrier(MockIbarrierRequest::NullaryPredicate cond) {
	  ibarrier_conditions.push_back(cond);
	}
      private:
	int rank, size;
	
	struct SendData {
	  const void* expected_buf;
	  int count;
	  MPI_Datatype type;
	  int dest;
	  int tag;
	};
	
	struct RecvData {
	  const void* ans_buf;
	  int count;
	  MPI_Datatype type;
	  int src;
	  int tag;
	};
	
	struct GatherData {
	  const void* send_buf;
	  int sendcount;
	  MPI_Datatype sendtype;
	  const void* ans_buf;
	  int recvcount;
	  MPI_Datatype recvtype;
	  int root;
	};
	
	struct GatherVData {
	  const void* send_buf;
	  int sendcount;
	  MPI_Datatype sendtype;
	  const void* ans_buf;
	  const int* recvcounts;
	  const std::vector<int> displs;
	  MPI_Datatype recvtype;
	  
	  int root;
	};
	
	mutable std::deque<SendData> send_data;
	mutable std::deque<SendData> ssend_data;
	mutable std::deque<RecvData> recv_data;
	mutable std::deque<GatherData> gather_data;
	mutable std::deque<GatherVData > gatherv_data;
	mutable std::deque<MockIbarrierRequest::NullaryPredicate> ibarrier_conditions;
	
	virtual void BcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IbcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; return nullptr; }
	virtual void AllreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IallreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; return nullptr; }
	virtual std::shared_ptr<Request> IreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; return nullptr; }
	virtual void ReduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; }
	
	virtual void GatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				void* recv, int recvcount, MPI_Datatype recvtype,
				int root) const {
	  CPPUNIT_ASSERT(gather_data.size());
	  
	  auto& gd = gather_data.front();
	  CPPUNIT_ASSERT_EQUAL(gd.sendcount, sendcount);
	  CPPUNIT_ASSERT_EQUAL(gd.sendtype, sendtype);
	  CPPUNIT_ASSERT_EQUAL(gd.root, root);
	  int elSize;
	  MPI_Type_size(gd.sendtype, &elSize);
	  
	  { // assert that send[:] == relevant bit of stored data
	    int size_bytes = sendcount * elSize;
	    const char* expected = static_cast<const char*>(gd.send_buf);
	    const char* actual = static_cast<const char*>(send);
	    for (auto i = 0; i < size_bytes; ++i)
	      CPPUNIT_ASSERT_EQUAL(expected[i], actual[i]);
	  }
	  
	  if (Rank() == root) {
	    CPPUNIT_ASSERT_EQUAL(gd.recvcount, recvcount);
	    CPPUNIT_ASSERT_EQUAL(gd.recvtype, recvtype);

	    std::memcpy(recv, gd.ans_buf, recvcount * elSize * Size());
	  }
	  gather_data.pop_front();
	}
	
	virtual void GathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype,
				 int root) const {
	  CPPUNIT_ASSERT(gatherv_data.size());

	  auto& gd = gatherv_data.front();
	  CPPUNIT_ASSERT_EQUAL(gd.sendcount, sendcount);
	  CPPUNIT_ASSERT_EQUAL(gd.sendtype, sendtype);
	  CPPUNIT_ASSERT_EQUAL(gd.root, root);
	  
	  int elSize;
	  MPI_Type_size(gd.sendtype, &elSize);
	  
	  { // assert that send[:] == relevant bit of stored data
	    int size_bytes = sendcount * elSize;
	    const char* expected = static_cast<const char*>(gd.send_buf);
	    const char* actual = static_cast<const char*>(sendbuf);
	    for (auto i = 0; i < size_bytes; ++i)
	      CPPUNIT_ASSERT_EQUAL(expected[i], actual[i]);
	  }
	  
	  if (Rank() == root) {
	    // The recv counts and displacements must match on root
	    for(auto i = 0; i < Size(); ++i) {
	      CPPUNIT_ASSERT_EQUAL(gd.recvcounts[i], recvcounts[i]);
	      CPPUNIT_ASSERT_EQUAL(gd.displs[i], displs[i]);
	    }
	    CPPUNIT_ASSERT_EQUAL(gd.recvtype, recvtype);
	    
	    int nElTotal = std::accumulate(recvcounts, recvcounts+Size(), 0);
	    std::memcpy(recvbuf, gd.ans_buf, nElTotal*elSize);
	  }
	  gatherv_data.pop_front();
	}
	
	virtual void AllgatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				   void* recv, int recvcount, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
	virtual void AllgathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				    void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
	virtual void AlltoallImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				  void* recv, int recvcount, MPI_Datatype recvtype) const { NOTIMPLEMENTED; }
	virtual void SendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			      int dest, int tag) const {
	  CPPUNIT_ASSERT(send_data.size());
	  
	  auto& sd = send_data.front();
	  CPPUNIT_ASSERT_EQUAL(sd.count, sendcount);
	  CPPUNIT_ASSERT_EQUAL(sd.type, sendtype);
	  CPPUNIT_ASSERT_EQUAL(sd.dest, dest);
	  CPPUNIT_ASSERT_EQUAL(sd.tag, tag);
	  send_data.pop_front();
	}
	virtual void SsendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			       int dest, int tag) const {
	  CPPUNIT_ASSERT(ssend_data.size());
	  
	  auto& sd = ssend_data.front();
	  CPPUNIT_ASSERT_EQUAL(sd.count, sendcount);
	  CPPUNIT_ASSERT_EQUAL(sd.type, sendtype);
	  CPPUNIT_ASSERT_EQUAL(sd.dest, dest);
	  CPPUNIT_ASSERT_EQUAL(sd.tag, tag);
	  ssend_data.pop_front();
	}
	
	virtual void RecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
			      int src, int tag, MPI_Status* stat) const {
	  CPPUNIT_ASSERT(recv_data.size());
	  
	  auto& rd = recv_data.front();
	  CPPUNIT_ASSERT_EQUAL(rd.count, recvcount);
	  CPPUNIT_ASSERT_EQUAL(rd.type, recvtype);
	  CPPUNIT_ASSERT_EQUAL(rd.src, src);
	  CPPUNIT_ASSERT_EQUAL(rd.tag, tag);
	  
	  int elSize;
	  MPI_Type_size(rd.type, &elSize);
	  
	  std::memcpy(recvbuf, rd.ans_buf, recvcount * elSize);
	  recv_data.pop_front();

	}
	
	virtual std::shared_ptr<Request> IsendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						   int dest, int tag) const {
	  SendImpl(sendbuf, sendcount, sendtype, dest, tag);
	  return std::make_shared<MockRequest>();
	}
	
	virtual std::shared_ptr<Request> IssendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						    int dest, int tag) const {
	  SsendImpl(sendbuf, sendcount, sendtype, dest, tag);
	  return std::make_shared<MockRequest>();
	}
	virtual std::shared_ptr<Request> IrecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
						   int source, int tag) const {
	  RecvImpl(recvbuf, recvcount, recvtype, source, tag, MPI_STATUS_IGNORE);
	  return std::make_shared<MockRequest>();
	}
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_MOCKCOMMUNICATOR_H
