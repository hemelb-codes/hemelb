// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_NET_MPICOMMUNICATOR_HPP
#define HEMELB_NET_MPICOMMUNICATOR_HPP

#include <numeric>
#include "net/MpiDataType.h"
#include "net/MpiConstness.h"

namespace hemelb
{
  namespace net
  {
    template<typename T>
    void MpiCommunicator::Broadcast(T& val, const int root) const
    {
      HEMELB_MPI_CALL(MPI_Bcast, (&val, 1, MpiDataType<T>(), root, *this));
    }
    template<typename T>
    void MpiCommunicator::Broadcast(std::vector<T>& vals, const int root) const
    {
      HEMELB_MPI_CALL(MPI_Bcast, (&vals[0], vals.size(), MpiDataType<T>(), root, *this));
    }

    template<typename T>
    T MpiCommunicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(MPI_Allreduce, (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(MPI_Allreduce,
                      (MpiConstCast(&vals[0]), &ans[0], vals.size(), MpiDataType<T>(), op, *this));
      return ans;
    }

    template<typename T>
    T MpiCommunicator::Reduce(const T& val, const MPI_Op& op, const int root) const
    {
      T ans;
      HEMELB_MPI_CALL(MPI_Reduce, (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, root, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Reduce(const std::vector<T>& vals, const MPI_Op& op,
                                           const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(vals.size());
        recvbuf = &ans[0];
      }

      HEMELB_MPI_CALL(MPI_Reduce,
                      (MpiConstCast(&vals[0]), recvbuf, vals.size(), MpiDataType<T>(), op, root, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const T& val, const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(Size());
        recvbuf = &ans[0];
      }
      HEMELB_MPI_CALL(MPI_Gather,
                      (MpiConstCast(&val), 1, MpiDataType<T>(), recvbuf, 1, MpiDataType<T>(), root, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const std::vector<T>& val, const int root) const
    {
      std::vector<T> ans;
      std::vector<int> offsets;
      auto const counts = Gather(static_cast<int>(val.size()), root);

      if (Rank() == root)
      {
        offsets.push_back(0);
        for(auto const &c: counts)
        {
          offsets.push_back(offsets.back() + c);
        }

        auto const N = std::accumulate(counts.begin(), counts.end(), 0);
        ans.reserve(N+1);
        ans.resize(N);
      }

      HEMELB_MPI_CALL(MPI_Gatherv,
                      (MpiConstCast(val.data()), val.size(), MpiDataType<T>(),
                       Rank() == root ? ans.data(): nullptr,
                       Rank() == root ? counts.data(): nullptr,
                       Rank() == root ? offsets.data(): nullptr,
                       MpiDataType<T>(), root, *this));
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

    template<typename T>
    std::vector<T> MpiCommunicator::AllGather(const T& val) const
    {
      std::vector<T> ans(Size());
      T* recvbuf = &ans[0];

      HEMELB_MPI_CALL(MPI_Allgather,
                      (MpiConstCast(&val), 1, MpiDataType<T>(), recvbuf, 1, MpiDataType<T>(), *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllNeighGather(const T& val) const
    {
      std::vector<T> ans(GetNeighborsCount());
      T* recvbuf = ans.data();

      HEMELB_MPI_CALL(MPI_Neighbor_allgather,
                      (MpiConstCast(&val), 1, MpiDataType<T>(), recvbuf, 1, MpiDataType<T>(), *this));

      return ans;
    }

    template<typename T>
    std::vector<std::vector<T>> MpiCommunicator::AllNeighGatherV(const std::vector<T>& val) const
    {
      int numProcs = GetNeighborsCount();
      std::vector<int> valSizes = AllNeighGather((int) val.size());
      std::vector<int> valDisplacements(numProcs + 1);

      int totalSize = std::accumulate(valSizes.begin(),
                                      valSizes.end(),
                                      0);

      valDisplacements[0] = 0;
      for (int j = 0; j < numProcs; ++j)
      {
        valDisplacements[j + 1] = valDisplacements[j] + valSizes[j];
      }

      std::vector<T> allVal(totalSize);
      HEMELB_MPI_CALL(MPI_Neighbor_allgatherv,
                      ( net::MpiConstCast(&val[0]), val.size(), net::MpiDataType<T>(), &allVal[0], net::MpiConstCast(&valSizes[0]), net::MpiConstCast(&valDisplacements[0]), net::MpiDataType<T>(), *this ));

      std::vector<std::vector<T>> ans(numProcs);
      for (decltype(numProcs) procIndex = 0; procIndex < numProcs; ++procIndex)
      {
        ans[procIndex].reserve(valDisplacements[procIndex + 1] - valDisplacements[procIndex]);
        for (auto indexAllCoords = valDisplacements[procIndex];
             indexAllCoords < valDisplacements[procIndex + 1]; ++indexAllCoords)
        {
          ans[procIndex].push_back(allVal[indexAllCoords]);
        }
      }

      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllToAll(const std::vector<T>& vals) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(MPI_Alltoall,
                      (MpiConstCast(&vals[0]), 1, MpiDataType<T>(), &ans[0], 1, MpiDataType<T>(), *this));
      return ans;
    }

    template<typename T>
    void MpiCommunicator::Send(const T& val, int dest, int tag) const
    {
      HEMELB_MPI_CALL(MPI_Send, (MpiConstCast(&val), 1, MpiDataType<T>(), dest, tag, *this));
    }
    template<typename T>
    void MpiCommunicator::Send(const std::vector<T>& vals, int dest, int tag) const
    {
      HEMELB_MPI_CALL(MPI_Send,
                      (MpiConstCast(&vals[0]), vals.size(), MpiDataType<T>(), dest, tag, *this));
    }

    template<typename T>
    void MpiCommunicator::Receive(T& val, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(MPI_Recv, (&val, 1, MpiDataType<T>(), src, tag, *this, stat));
    }
    template<typename T>
    void MpiCommunicator::Receive(std::vector<T>& vals, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(MPI_Recv, (&vals, vals.size(), MpiDataType<T>(), src, tag, *this, stat));
    }
  }
}

#endif // HEMELB_NET_MPICOMMUNICATOR_HPP
