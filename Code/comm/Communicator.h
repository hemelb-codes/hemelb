// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_COMMUNICATOR_H
#define HEMELB_COMM_COMMUNICATOR_H

#include <vector>
#include <memory>
#include <mpi.h>

namespace hemelb
{
  namespace comm
  {
    class Group;
    class Request;
    class RequestList;
    class MpiFile;

    // Base class for communicators (MPI, null, and mock)
    class Communicator
    {
      public:
        typedef std::shared_ptr<Communicator> Ptr;
        typedef std::shared_ptr<const Communicator> ConstPtr;
        /**
         * Class has virtual methods so should have virtual d'tor.
         */
        virtual ~Communicator() {};

        /**
         * Returns the local rank on the communicator
         * @return
         */
        virtual int Rank() const = 0;
      
        inline bool OnIORank() const
        {
	  return Rank() == GetIORank();
	}
      
        inline int GetIORank() const
        {
	  return 0;
	}

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        virtual int Size() const = 0;

        /**
         * Abort - should try to bring down all tasks, but no guarantees
         * @param errCode
         */
        virtual void Abort(int errCode) const = 0;

        /**
         * Duplicate the communicator - see MPI_COMM_DUP
         * @return
         */
        virtual Ptr Duplicate() const = 0;

        virtual std::shared_ptr<Group> GetGroup() const = 0;
        virtual Ptr Create(std::shared_ptr<const Group> grp) const = 0;
        /**
         * Opens a file with MPI_File_open. A collective operation
         * @param filename
         * @param mode
         * @param info
         * @return
         */
        virtual MpiFile OpenFile(const std::string& filename, int mode,
				 const MPI_Info info = MPI_INFO_NULL) const = 0;

        virtual std::shared_ptr<RequestList> MakeRequestList() const = 0;
        virtual void Barrier() const = 0;
        virtual std::shared_ptr<Request> Ibarrier() const = 0;

        virtual bool Iprobe(int source, int tag, MPI_Status* stat=MPI_STATUS_IGNORE) const = 0;

        template <typename T>
        void Broadcast(T& val, const int root) const;
        template <typename T>
        void Broadcast(std::vector<T>& vals, const int root) const;

        template <typename T>
        std::shared_ptr<Request> Ibcast(T& val, const int root) const;
        template <typename T>
        std::shared_ptr<Request> Ibcast(std::vector<T>& vals, const int root) const;

        template <typename T>
        T AllReduce(const T& val, const MPI_Op& op) const;
        template <typename T>
        std::vector<T> AllReduce(const std::vector<T>& vals, const MPI_Op& op) const;

        template <typename T>
        std::shared_ptr<Request> Iallreduce(const T& val, const MPI_Op& op, T& out) const;

        template <typename T>
        std::shared_ptr<Request> Ireduce(const T& val, const MPI_Op& op, const int root, T& out) const;

        template <typename T>
        T Reduce(const T& val, const MPI_Op& op, const int root) const;
        template <typename T>
        std::vector<T> Reduce(const std::vector<T>& vals, const MPI_Op& op, const int root) const;

        template <typename T>
        std::vector<T> Gather(const T& val, const int root) const;
        
        template <typename T>
        std::vector<T> GatherV(const std::vector<T> senddata, const std::vector<int> recvcounts,
			       const int root) const;
        template <typename T>
        std::vector<T> AllGatherV(const std::vector<T> senddata, const std::vector<int> recvcounts) const;
	      
        template <typename T>
        std::vector<T> AllGather(const T& val) const;

        template <typename T>
        T Scatter(const std::vector<T>& vals, const int root) const;
        template <typename T>
	std::vector<T> Scatter(const std::vector<T>& vals, const size_t n, const int root) const;

        template <typename T>
        std::vector<T> AllToAll(const std::vector<T>& vals) const;

        template <typename T>
        void Send(const std::vector<T>& val, int dest, int tag=0) const;
        template <typename T>
        void Send(const T& val, int dest, int tag=0) const;
        template <typename T>
        void Send(const T* valPtr, int count, int dest, int tag) const;

        template <typename T>
        void Recv(std::vector<T>& val, int src, int tag=0, MPI_Status* stat=MPI_STATUS_IGNORE) const;
        template <typename T>
        void Recv(T& val, int src, int tag=0, MPI_Status* stat=MPI_STATUS_IGNORE) const;
        template <typename T>
        void Recv(T* val, int count, int src, int tag=0, MPI_Status* stat=MPI_STATUS_IGNORE) const;

        template <typename T>
        std::shared_ptr<Request> Isend(const T& val, int dest, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Isend(const std::vector<T>& vals, int dest, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Isend(const T* valPtr, int count, int dest, int tag=0) const;

        template <typename T>
        std::shared_ptr<Request> Issend(const T& val, int dest, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Issend(const std::vector<T>& vals, int dest, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Issend(const T* valPtr, int count, int dest, int tag=0) const;

        template <typename T>
        std::shared_ptr<Request> Irecv(T& val, int source, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Irecv(std::vector<T>& vals, int source, int tag=0) const;
        template <typename T>
        std::shared_ptr<Request> Irecv(T* valPtr, int count, int source, int tag=0) const;


    protected:
      virtual void BcastImpl(void* buf, int count, MPI_Datatype dt, int root) const = 0;
      virtual std::shared_ptr<Request> IbcastImpl(void* buf, int count, MPI_Datatype dt, int root) const = 0;
      virtual void AllreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const = 0;
      virtual std::shared_ptr<Request> IallreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op) const = 0;
      virtual std::shared_ptr<Request> IreduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const = 0;
      virtual void ReduceImpl(const void* send, void* ans, int count, MPI_Datatype dt, MPI_Op op, int root) const = 0;
      virtual void GatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
			      void* recv, int recvcount, MPI_Datatype recvtype,
			      int root) const = 0;
      virtual void GathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			       void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype,
			       int root) const = 0;
      virtual void AllgatherImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				 void* recv, int recvcount, MPI_Datatype recvtype) const = 0;
      virtual void AllgathervImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
				  void *recvbuf, const int* recvcounts, const int* displs, MPI_Datatype recvtype) const = 0;
      virtual void ScatterImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			       void* recvbuf, int recvcount, MPI_Datatype recvtype, int root) const = 0;
      virtual void AlltoallImpl(const void* send, int sendcount, MPI_Datatype sendtype,
				void* recv, int recvcount, MPI_Datatype recvtype) const = 0;
      virtual void SendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			    int dest, int tag) const = 0;
      virtual void SsendImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			     int dest, int tag) const = 0;
      virtual void RecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
			    int src, int tag, MPI_Status* stat) const = 0;
      virtual std::shared_ptr<Request> IsendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
				 int dest, int tag) const = 0;
      virtual std::shared_ptr<Request> IssendImpl(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
				  int dest, int tag) const = 0;
      virtual std::shared_ptr<Request> IrecvImpl(void* recvbuf, int recvcount, MPI_Datatype recvtype,
				 int source, int tag) const = 0;
    };

    // bool operator==(const Communicator& comm1, const Communicator& comm2);
    // bool operator!=(const Communicator& comm1, const Communicator& comm2);

    template <class T>
    struct is_std_vector : public std::false_type {};
    template <class T>
    struct is_std_vector<std::vector<T> > : public std::true_type {};
  }
}

#include "comm/Communicator.hpp"

#endif /* HEMELB_COMM_MPICOMMUNICATOR_H */
