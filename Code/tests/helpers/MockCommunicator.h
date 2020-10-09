
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_MOCKCOMMUNICATOR_H
#define HEMELB_TESTS_HELPERS_MOCKCOMMUNICATOR_H

#include <deque>
#include <numeric>
#include <cstring>

#include <catch2/catch.hpp>

#include "comm/Communicator.h"
#include "comm/MpiFile.h"
#include "comm/Request.h"

namespace hemelb
{
  namespace tests
  {
    namespace helpers
    {
      using namespace comm;

      class MockRequest : public comm::Request
      {
      public:
	MockRequest();
	virtual ~MockRequest() = default;
        virtual void Wait();
        virtual bool Test();

      private:
	int testCount;
      };

      class MockIbarrierRequest : public MockRequest
      {
      public:
	typedef std::function<bool()> NullaryPredicate;
	MockIbarrierRequest(NullaryPredicate fn);
	virtual ~MockIbarrierRequest() = default;
	virtual void Wait();
	virtual bool Test();
      private:
	NullaryPredicate done;
      };
      
      class MockRequestList : public comm::RequestList
      {
	std::vector<Request::Ptr> reqs;
      public:
	virtual ~MockRequestList() = default;

	virtual size_t size() const;
	virtual void resize(size_t i);
	virtual void push_back(Request::Ptr rp);
	
	virtual void clear();
	
	virtual void set(size_t i, Request::Ptr rp);
	
	virtual void WaitAll();
	virtual bool TestAll();
      };
            
      class MockCommunicator : public comm::Communicator
      {
      public:
	/***
	 * Constructor for a dummy communicator
	 * Can be useful for testing but can't actually be used
	 * @param rank
	 * @param size
	 */
	MockCommunicator(int rank_, int size_);

	virtual int Rank() const;
	virtual int Size() const;
	
	virtual void Abort(int errCode) const;
	
        virtual Ptr Duplicate() const;

        virtual std::shared_ptr<Group> GetGroup() const;
        virtual Ptr Create(std::shared_ptr<const Group> grp) const;
        virtual MpiFile OpenFile(const std::string& filename, int mode,
				 const MPI_Info info = MPI_INFO_NULL) const;
      
        virtual std::shared_ptr<RequestList> MakeRequestList() const;
        virtual void Barrier() const;
	
        virtual std::shared_ptr<Request> Ibarrier() const;

        virtual bool Iprobe(int source, int tag, MPI_Status* stat=MPI_STATUS_IGNORE) const;

	void ExpectationsAllCompleted();
	
	template<class T>
	void RequireGather(const std::vector<T>& data, int root)
	{
	  INFO("Number of elements must be a multiple of comm size");
	  REQUIRE(0 == int(data.size()) % Size());
	  int sendcount = data.size() / Size();
	  gather_data.push_back({data.data() + sendcount*Rank(), sendcount, comm::MpiDataType<T>(),
		data.data(), sendcount, comm::MpiDataType<T>(),
		root});
	}
	
	template<class T>
	void RequireGatherV(const std::vector<T>& data, const std::vector<int>& recvcounts, int root)
	{
	  REQUIRE(Size() == int(recvcounts.size()));
	  std::vector<int> displs(Size());
	  int total = 0;
	  for(size_t i = 0; i < Size(); ++i) {
	    displs[i] = total;
	    total += recvcounts[i];
	  }
	  auto my_send_data = data.data() + displs[Rank()];
	  INFO("Size of data does not match sum of recvcounts");
	  REQUIRE(size_t(total) == data.size());
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

	  REQUIRE(dest >= 0);
	  REQUIRE(dest < Size());
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
	  REQUIRE(src >= 0);
	  REQUIRE(src < Size());
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

	void RequireIbarrier(MockIbarrierRequest::NullaryPredicate cond);
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
	
	virtual void BcastImpl(void* buf, int count, MPI_Datatype dt, int root) const;
	virtual std::shared_ptr<Request> IbcastImpl(void* buf, int count, MPI_Datatype dt, int root) const;
	virtual void AllreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const;
	virtual std::shared_ptr<Request> IallreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const;
	virtual std::shared_ptr<Request> IreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const;
	virtual void ReduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const;
	
	virtual void GatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				void* recv, int recvcount, MPI_Datatype recvtype,
				int root) const;
	
	virtual void GathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype,
				 int root) const;
	
	virtual void AllgatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				   void* recv, int recvcount, MPI_Datatype recvtype) const;
	virtual void AllgathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				    void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype) const;
	virtual void ScatterImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				 void* recvbuf, int recvcount, MPI_Datatype recvtype, int root) const;
	virtual void AlltoallImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				  void* recv, int recvcount, MPI_Datatype recvtype) const;
	virtual void SendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			      int dest, int tag) const;
	virtual void SsendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			       int dest, int tag) const;
	
	virtual void RecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
			      int src, int tag, MPI_Status* stat) const;
	
	virtual std::shared_ptr<Request> IsendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						   int dest, int tag) const;
	
	virtual std::shared_ptr<Request> IssendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						    int dest, int tag) const;
	virtual std::shared_ptr<Request> IrecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
						   int source, int tag) const;
      };

    }
  }
}

#endif // HEMELB_TESTS_HELPERS_MOCKCOMMUNICATOR_H
