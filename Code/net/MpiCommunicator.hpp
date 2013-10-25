//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_NET_MPICOMMUNICATOR_HPP
#define HEMELB_NET_MPICOMMUNICATOR_HPP

#include "net/MpiDataType.h"

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
    T MpiCommunicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (const_cast<T*>(&val), &ans, 1, MpiDataType<T>(), op, *this)
      );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (const_cast<T*>(&vals[0]), &ans[0], vals.size(), MpiDataType<T>(), op, *this)
      );
      return ans;
    }

  }
}

#endif // HEMELB_NET_MPICOMMUNICATOR_HPP
