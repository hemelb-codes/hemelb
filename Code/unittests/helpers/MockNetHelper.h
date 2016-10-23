
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H
#define HEMELB_UNITTESTS_HELPERS_MOCKNETHELPER_H

#include "unittests/net/NetMock.h"
#include <numeric>

namespace hemelb
{
  namespace unittests
  {
    namespace helpers
    {
      using namespace comm;

      class MockMpiCommunicator : public comm::Communicator
      {
      public:
	/***
	 * Constructor for a dummy communicator
	 * Can be useful for testing but can't actually be used
	 * @param rank
	 * @param size
	 */
	MockMpiCommunicator(int rank_, int size_) :
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
	
        virtual Ptr Duplicate() const { NOTIMPLEMENTED; }

        virtual std::shared_ptr<Group> GetGroup() const { NOTIMPLEMENTED; }
        virtual Ptr Create(std::shared_ptr<const Group> grp) const { NOTIMPLEMENTED; }
        virtual std::shared_ptr<MpiFile> OpenFile(const std::string& filename, int mode,
						  const MPI_Info info = MPI_INFO_NULL) const { NOTIMPLEMENTED; }
      
        virtual std::shared_ptr<RequestList> MakeRequestList() const { NOTIMPLEMENTED; }
        virtual void Barrier() const { NOTIMPLEMENTED; }
        virtual std::shared_ptr<Request> Ibarrier() const { NOTIMPLEMENTED; }

        virtual bool Iprobe(int source, int tag, MPI_Status* stat=MPI_STATUS_IGNORE) const { NOTIMPLEMENTED; }
	template<class T>
	void AddGatherResult(const std::vector<T>& data)
	{
	  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of elements must be a multiple of comm size",
				       0, int(data.size()) % Size());
	  gather_data.emplace_back(data.size() / Size(), sizeof(T), data.data());
	}
	
	template<class T>
	void AddGatherVResult(const std::vector<T>& data)
	{
	  gatherv_data.emplace_back(data.size(), sizeof(T), data.data());
	}
      private:
	int rank, size;
	typedef std::tuple<int, int, const void*> gather_info;
	mutable std::deque<gather_info > gather_data;
	mutable std::deque<gather_info > gatherv_data;
	
	virtual void BcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IbcastImpl(void* buf, int count, MPI_Datatype dt, int root) const { NOTIMPLEMENTED; }
	virtual void AllreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IallreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; }
	virtual void ReduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const { NOTIMPLEMENTED; }
	
	virtual void GatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				void* recv, int recvcount, MPI_Datatype recvtype,
				int root) const {
	  CPPUNIT_ASSERT(gather_data.size());
	  int nEl, elSize;
	  const void* ptr;
	  std::tie(nEl, elSize, ptr) = gather_data.front();
	  CPPUNIT_ASSERT_EQUAL(nEl, sendcount);
	  
	  { // assert that send[:] == stored_data[sendcount*rank:sendcount*(rank+1)]
	    int size_bytes = nEl*elSize;
	    const char* expected = static_cast<const char*>(ptr) + Rank()*size_bytes;
	    const char* actual = static_cast<const char*>(send);
	    for (auto i = 0; i < size_bytes; ++i)
	      CPPUNIT_ASSERT_EQUAL(expected[i], actual[i]);
	  }
	  
	  if (Rank() == root) {
	    std::memcpy(recv, ptr, nEl*elSize*Size());
	  }
	  gather_data.pop_front();
	}
	
	virtual void GathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype,
				 int root) const {
	  CPPUNIT_ASSERT(gatherv_data.size());
	  int nElTotal, elSize;
	  const void* ptr;
	  std::tie(nElTotal, elSize, ptr) = gatherv_data.front();
	  if (Rank() == root)
	    CPPUNIT_ASSERT_EQUAL_MESSAGE("This rank's count must match",
					 recvcounts[Rank()], sendcount);

	  if (root == Rank()) {
	    // sum(recvcounts) == nElTotal (only significant on root of collective)
	    int sum_recvcounts = std::accumulate(recvcounts, recvcounts + Size(), 0);
	    CPPUNIT_ASSERT_EQUAL(nElTotal, sum_recvcounts);
	  }
	  // assert that send[:] == stored_data[sendcount*rank:sendcount*(rank+1)]
	  if (Rank() == root) {
	    std::memcpy(recvbuf, ptr, nElTotal*elSize);
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
			      int dest, int tag) const { NOTIMPLEMENTED; }
	virtual void SsendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			       int dest, int tag) const { NOTIMPLEMENTED; }
	virtual void RecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
			      int src, int tag, MPI_Status* stat) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IsendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						   int dest, int tag) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IssendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
						    int dest, int tag) const { NOTIMPLEMENTED; }
	virtual std::shared_ptr<Request> IrecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
						   int source, int tag) const { NOTIMPLEMENTED; }
      };

      class MockNetHelper
      {
      protected:
	MockNetHelper() :
	  communicatorMock(NULL), netMock(NULL)
	{
	}
	void setUp(const proc_t core_count, const proc_t current_core)
	{
	  communicatorMock = std::make_shared<MockMpiCommunicator>(current_core, core_count);
	  netMock = new net::NetMock(communicatorMock);
	}
	void tearDown()
	{
	  delete netMock;
	}

	comm::Communicator::Ptr communicatorMock;
	net::NetMock *netMock;
      };

    }
  }
}

#endif // HEMELB_UNITTESTS_HELPERS_RANDOMSOURCE_H
