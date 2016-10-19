
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_COMM_COMMUNICATOR_HPP
#define HEMELB_COMM_COMMUNICATOR_HPP

#include "comm/MpiDataType.h"

namespace hemelb
{
  namespace comm
  {
    template<typename T>
    void Communicator::Broadcast(T& val, const int root) const
    {
      BcastImpl(&val, 1, MpiDataType<T>(), root);
    }
    template<typename T>
    void Communicator::Broadcast(std::vector<T>& vals, const int root) const
    {
      BcastImpl(vals.data(), vals.size(), MpiDataType<T>(), root);
    }

    template <typename T>
    std::shared_ptr<Request> Communicator::Ibcast(T& val, const int root) const
    {
      return IbcastImpl(&val, 1, MpiDataType<T>(), root);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Ibcast(std::vector<T>& vals, const int root) const
    {
      return IbcastImpl(vals.data(), vals.size(), MpiDataType<T>(), root);
    }

    template<typename T>
    T Communicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      AllreduceImpl(&val, &ans, 1, MpiDataType<T>(), op);
      return ans;
    }

    template<typename T>
    std::vector<T> Communicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      AllreduceImpl(vals.data(), ans.data(), vals.size(), MpiDataType<T>(), op);
      return ans;
    }

    template <typename T>
    std::shared_ptr<Request> Communicator::Iallreduce(const T& val, const MPI_Op& op, T& out) const
    {
      return IallreduceImpl(&val, &out, 1, MpiDataType<T>(), op);
    }

    template <typename T>
    std::shared_ptr<Request> Communicator::Ireduce(const T& val, const MPI_Op& op, const int root, T& out) const
    {
      return IreduceImpl(&val, &out, 1, MpiDataType<T>(), op, root);
    }

    template<typename T>
    T Communicator::Reduce(const T& val, const MPI_Op& op, const int root) const
    {
      T ans;
      ReduceImpl(&val, &ans, 1, MpiDataType<T>(), op, root);
      return ans;
    }

    template<typename T>
    std::vector<T> Communicator::Reduce(const std::vector<T>& vals, const MPI_Op& op,
					const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(vals.size());
        recvbuf = ans.data();
      }

      ReduceImpl(vals.data(), recvbuf, vals.size(), MpiDataType<T>(), op, root);
      return ans;
    }

    template<typename T>
    std::vector<T> Communicator::Gather(const T& val, const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(Size());
        recvbuf = ans.data();
      }
      
      GatherImpl(&val, 1, MpiDataType<T>(),
		 recvbuf, 1, MpiDataType<T>(),
		 root);
      return ans;
    }

    template <typename T>
    std::vector<T> Communicator::GatherV(const std::vector<T> senddata,
					    const std::vector<int> recvcounts,
					    const int root) const
    {
      const int np = Size();
      const int sendcount = senddata.size();
      std::vector<int> displs;
      std::vector<T> ans;
      T* recvbuf = nullptr;
      if (Rank() == root)
      {
	// Compute the displacements from the counts
	displs.resize(np);
	int total = 0;
	for(size_t i = 0; i < np; ++i) {
	  displs[i] = total;
	  total += recvcounts[i];
	}
	// set up recv buffer
	ans.resize(total);
	recvbuf = ans.data();
      }

      GathervImpl(senddata.data(), sendcount, MpiDataType<T>(),
		  recvbuf, recvcounts.data(), displs.data(), MpiDataType<T>(),
		  root);
      return ans;
    }
    
    template<typename T>
    std::vector<T> Communicator::AllGather(const T& val) const
    {
      std::vector<T> ans(Size());
      T* recvbuf =  ans.data();
      
      AllgatherImpl(&val, 1, MpiDataType<T>(),
		    recvbuf, 1, MpiDataType<T>());
      return ans;
    }

    template <typename T>
    std::vector<T> Communicator::AllToAll(const std::vector<T>& vals) const
    {
      std::vector<T> ans(vals.size());
      AlltoallImpl(vals.data(), 1, MpiDataType<T>(),
		   ans.data(), 1, MpiDataType<T>());
      return ans;
    }

    // Send implementations
    template <typename T>
    void Communicator::Send(const T* valPtr, int count, int dest, int tag) const
    {
      SendImpl(valPtr, count, MpiDataType<T>(), dest, tag);
    }
    template <typename T>
    void Communicator::Send(const T& val, int dest, int tag) const
    {
      Send(&val, 1, dest, tag);
    }
    template <typename T>
    void Communicator::Send(const std::vector<T>& vals, int dest, int tag) const
    {
      Send(vals.data(), vals.size(), dest, tag);
    }

    // Recv implementations
    template <typename T>
    void Communicator::Recv(T* valPtr, int count, int src, int tag, MPI_Status* stat) const
    {
      RecvImpl(valPtr, count, MpiDataType<T>(), src, tag, stat);
    }
    template <typename T>
    void Communicator::Recv(T& val, int src, int tag, MPI_Status* stat) const
    {
      Recv(&val, 1, src, tag, stat);
    }
    template <typename T>
    void Communicator::Recv(std::vector<T>& vals, int src, int tag, MPI_Status* stat) const
    {
      Recv(vals.data(), vals.size(), src, tag, stat);
    }

    // Isend implementations
    template <typename T>
    std::shared_ptr<Request> Communicator::Isend(const T* valPtr, int count, int dest, int tag) const
    {
      return IsendImpl(valPtr, count, MpiDataType<T>(), dest, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Isend(const T& val, int dest, int tag) const
    {
      return Isend(&val, 1, dest, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Isend(const std::vector<T>& vals, int dest, int tag) const
    {
      return Isend(vals.data(), vals.size(), dest, tag);
    }

    // Issend implementations
    template <typename T>
    std::shared_ptr<Request> Communicator::Issend(const T* valPtr, int count, int dest, int tag) const
    {
     return IssendImpl(valPtr, count, MpiDataType<T>(), dest, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Issend(const T& val, int dest, int tag) const
    {
      return Issend(&val, 1, dest, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Issend(const std::vector<T>& vals, int dest, int tag) const
    {
      return Issend(vals.data(), vals.size(), dest, tag);
    }

    // Irecv implementations
    template <typename T>
    std::shared_ptr<Request> Communicator::Irecv(T* valPtr, int count, int source, int tag) const
    {
      return IrecvImpl(valPtr, count, MpiDataType<T>(), source, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Irecv(T& val, int source, int tag) const
    {
      return Irecv(&val, 1, source, tag);
    }
    template <typename T>
    std::shared_ptr<Request> Communicator::Irecv(std::vector<T>& vals, int source, int tag) const
    {
      return Irecv(&vals[0], vals.size(), source, tag);
    }

  }
}

#endif // HEMELB_COMM_COMMUNICATOR_HPP
