// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COMM_MPICOMMUNICATOR_H
#define HEMELB_COMM_MPICOMMUNICATOR_H

#include "comm/Communicator.h"

namespace hemelb
{
  namespace comm
  {
    
    // MPI communicator
    class MpiCommunicator : public Communicator, public std::enable_shared_from_this<MpiCommunicator>
    {
      public:
        /**
	 * Constructor for an uninitialised communicator, equivalent to
	 * MPI_COMM_NULL
	 * @param communicator
	 */
        MpiCommunicator();
        /**
	 * Constructor to get data needed from an MPI communicator
	 * @param communicator
	 */
        MpiCommunicator(MPI_Comm communicator, bool willOwn);

	operator MPI_Comm() const;
        /**
         * Class has virtual methods so should have virtual d'tor.
         */
        virtual ~MpiCommunicator();

        /**
         * Returns the local rank on the communicator
         * @return
         */
        virtual int Rank() const;

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        virtual int Size() const;

        /**
         * Abort - should try to bring down all tasks, but no guarantees
         * @param errCode
         */
        virtual void Abort(int errCode) const;

        /**
         * Duplicate the communicator - see MPI_COMM_DUP
         * @return
         */
      virtual Communicator::Ptr Duplicate() const;
	
      virtual std::shared_ptr<Group> GetGroup() const;
      virtual Communicator::Ptr Create(std::shared_ptr<const Group> grp) const;

      virtual std::shared_ptr<MpiFile> OpenFile(const std::string& filename, int mode,
				  const MPI_Info info = MPI_INFO_NULL) const;

        virtual void Barrier() const;
        virtual std::shared_ptr<Request> Ibarrier() const;

        virtual bool Iprobe(int source, int tag, MPI_Status* stat=MPI_STATUS_IGNORE) const;

    private:
      
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
      
      std::shared_ptr<MPI_Comm> commPtr;
    };

    // bool operator==(const Communicator& comm1, const Communicator& comm2);
    // bool operator!=(const Communicator& comm1, const Communicator& comm2);

  }
}


#endif /* HEMELB_COMM_MPICOMMUNICATOR_H */
