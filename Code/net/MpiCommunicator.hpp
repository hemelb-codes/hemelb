
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_NET_MPICOMMUNICATOR_HPP
#define HEMELB_NET_MPICOMMUNICATOR_HPP

#include "net/MpiDataType.h"
#include "net/MpiConstness.h"
#include "net/MpiRequest.h"

namespace hemelb
{
  namespace net
  {
    template<typename T>
    void MpiCommunicator::Broadcast(T& val, const int root) const
    {
      HEMELB_MPI_CALL(
          MPI_Bcast,
          (&val, 1, MpiDataType<T>(), root, *this)
      );
    }
    template<typename T>
    void MpiCommunicator::Broadcast(std::vector<T>& vals, const int root) const
    {
      HEMELB_MPI_CALL(
          MPI_Bcast,
          (&vals[0], vals.size(), MpiDataType<T>(), root, *this)
      );
    }

    template <typename T>
    MpiRequest MpiCommunicator::Ibcast(T& val, const int root) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Ibcast,
          (&val, 1, MpiDataType<T>(), root, *commPtr, &req)
      );
      return MpiRequest(req);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Ibcast(std::vector<T>& vals, const int root) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Ibcast,
          (&vals[0], vals.size(), MpiDataType<T>(), root, *commPtr, &req)
      );
      return MpiRequest(req);
    }

    template<typename T>
    T MpiCommunicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (MpiConstCast(&vals[0]), &ans[0], vals.size(), MpiDataType<T>(), op, *this)
      );
      return ans;
    }

    template <typename T>
    MpiRequest MpiCommunicator::Iallreduce(const T& val, const MPI_Op& op, T& out) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Iallreduce,
          (&val, &out, 1, MpiDataType<T>(), op, *commPtr, &req)
      );
      return MpiRequest(req);
    }

    template <typename T>
    MpiRequest MpiCommunicator::Ireduce(const T& val, const MPI_Op& op, const int root, T& out) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Ireduce,
          (&val, &out, 1, MpiDataType<T>(), op, root, *commPtr, &req)
      );
      return MpiRequest(req);
    }

    template<typename T>
    T MpiCommunicator::Reduce(const T& val, const MPI_Op& op, const int root) const
    {
      T ans;
      HEMELB_MPI_CALL(
          MPI_Reduce,
          (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, root, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Reduce(const std::vector<T>& vals, const MPI_Op& op,
                                           const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = NULL;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(vals.size());
        recvbuf = &ans[0];
      }

      HEMELB_MPI_CALL(
          MPI_Reduce,
          (MpiConstCast(&vals[0]), recvbuf, vals.size(), MpiDataType<T>(), op, root, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const T& val, const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = NULL;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(Size());
        recvbuf = &ans[0];
      }
      HEMELB_MPI_CALL(
          MPI_Gather,
          (MpiConstCast(&val), 1, MpiDataType<T>(),
              recvbuf, 1, MpiDataType<T>(),
              root, *this)
      );
      return ans;
    }

    template <typename T>
    T MpiCommunicator::Scatter(const std::vector<T>& vals, const int root) const {
      T ans;
      const T* ptr = (Rank() == root) ? vals.data() : nullptr;
      HEMELB_MPI_CALL(
		      MPI_Scatter,
		      (MpiConstCast(ptr), 1, MpiDataType<T>(),
		       &ans, 1, MpiDataType<T>(),
		       root, *this)
		      );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::Scatter(const std::vector<T>& vals, const size_t n, const int root) const {
      std::vector<T> ans(n);
      const T* ptr = (Rank() == root) ? vals.data() : nullptr;
      HEMELB_MPI_CALL(
		      MPI_Scatter,
		      (MpiConstCast(ptr), n, MpiDataType<T>(),
		       ans.data(), n, MpiDataType<T>(),
		       root, *this)
		      );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::GatherV(const std::vector<T> senddata,
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

      HEMELB_MPI_CALL(
	  MPI_Gatherv,
	  (MpiConstCast(senddata.data()), sendcount, MpiDataType<T>(),
	   recvbuf, recvcounts.data(), displs.data(), MpiDataType<T>(),
	   root, *this)
      );
      return ans;
    }
    
    template<typename T>
    std::vector<T> MpiCommunicator::AllGather(const T& val) const
    {
      std::vector<T> ans(Size());
      T* recvbuf =  &ans[0];

      HEMELB_MPI_CALL(
          MPI_Allgather,
          (MpiConstCast(&val), 1, MpiDataType<T>(),
              recvbuf, 1, MpiDataType<T>(),
              *this)
          );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::AllToAll(const std::vector<T>& vals) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(
          MPI_Alltoall,
          (MpiConstCast(&vals[0]), 1, MpiDataType<T>(),
           &ans[0], 1, MpiDataType<T>(),
           *this)
      );
      return ans;
    }

    // Send implementations
    template <typename T>
    void MpiCommunicator::Send(const T* valPtr, int count, int dest, int tag) const
    {
      HEMELB_MPI_CALL(
          MPI_Send,
          (valPtr, count, MpiDataType<T>(), dest, tag, *this)
      );
    }
    template <typename T>
    void MpiCommunicator::Send(const T& val, int dest, int tag) const
    {
      Send(&val, 1, dest, tag);
    }
    template <typename T>
    void MpiCommunicator::Send(const std::vector<T>& vals, int dest, int tag) const
    {
      Send(&vals[0], vals.size(), dest, tag);
    }

    // Recv implementations
    template <typename T>
    void MpiCommunicator::Recv(T* valPtr, int count, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(
          MPI_Recv,
          (valPtr, count, MpiDataType<T>(), src, tag, *this, stat)
      );
    }
    template <typename T>
    void MpiCommunicator::Recv(T& val, int src, int tag, MPI_Status* stat) const
    {
      Recv(&val, 1, src, tag);
    }
    template <typename T>
    void MpiCommunicator::Recv(std::vector<T>& vals, int src, int tag, MPI_Status* stat) const
    {
      Recv(&vals[0], vals.size(), src, tag);
    }

    // Isend implementations
    template <typename T>
    MpiRequest MpiCommunicator::Isend(const T* valPtr, int count, int dest, int tag) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Isend,
          (const_cast<T*>(valPtr), count, MpiDataType<T>(), dest, tag, *this, &req)
      );
      return MpiRequest(req);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Isend(const T& val, int dest, int tag) const
    {
      return Isend(&val, 1, dest, tag);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Isend(const std::vector<T>& vals, int dest, int tag) const
    {
      return Isend(&vals[0], vals.size(), dest, tag);
    }

    // Issend implementations
    template <typename T>
    MpiRequest MpiCommunicator::Issend(const T* valPtr, int count, int dest, int tag) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Issend,
          (const_cast<T*>(valPtr), count, MpiDataType<T>(), dest, tag, *this, &req)
      );
      return MpiRequest(req);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Issend(const T& val, int dest, int tag) const
    {
      return Issend(&val, 1, dest, tag);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Issend(const std::vector<T>& vals, int dest, int tag) const
    {
      return Issend(&vals[0], vals.size(), dest, tag);
    }

    // Irecv implementations
    template <typename T>
    MpiRequest MpiCommunicator::Irecv(T* valPtr, int count, int source, int tag) const
    {
      MPI_Request req;
      HEMELB_MPI_CALL(
          MPI_Irecv,
          (valPtr, count, MpiDataType<T>(), source, tag, *this, &req)
      );
      return MpiRequest(req);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Irecv(T& val, int source, int tag) const
    {
      return Irecv(&val, 1, source, tag);
    }
    template <typename T>
    MpiRequest MpiCommunicator::Irecv(std::vector<T>& vals, int source, int tag) const
    {
      return Irecv(&vals[0], vals.size(), source, tag);
    }

  }
}

#endif // HEMELB_NET_MPICOMMUNICATOR_HPP
