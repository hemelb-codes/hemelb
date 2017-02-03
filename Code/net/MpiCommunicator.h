// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_NET_MPICOMMUNICATOR_H
#define HEMELB_NET_MPICOMMUNICATOR_H

//#include "units.h"
//#include "net/mpi.h"
#include <vector>
#include "net/MpiError.h"
#include <memory>

namespace hemelb
{
  namespace net
  {
    class MpiGroup;

    class MpiCommunicator
    {
      public:
        static MpiCommunicator World();

        /**
         * Constructor for an uninitialised communicator, equivalent to
         * MPI_COMM_NULL
         * @param communicator
         */
        MpiCommunicator();

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
         * Creates a new communicator - see MPI_COMM_GROUP
         * @param Group which is a subset of the group of this communicator.
         * @return New communicator.
         */
        MpiCommunicator Create(const MpiGroup& grp) const;

        /**
         * Allow implicit casts to MPI_Comm
         * @return The underlying MPI communicator.
         */
        operator MPI_Comm() const
        {
          return *commPtr;
        }
        /**
         * Is this communicator valid? I.e. not equal to MPI_COMM_NULL.
         */
        operator bool() const
        {
          return (bool) commPtr;
        }
        /**
         * Returns the MPI group being used.
         * @return
         */
        MpiGroup Group() const;

        /**
         * Abort - should try to bring down all tasks, but no guarantees
         * @param errCode
         */
        void Abort(int errCode) const;

        /**
         * Duplicate the communicator - see MPI_COMM_DUP
         * @return
         */
        MpiCommunicator Duplicate() const;

        template<typename T>
        void Broadcast(T& val, const int root) const;
        template<typename T>
        void Broadcast(std::vector<T>& vals, const int root) const;

        template<typename T>
        T AllReduce(const T& val, const MPI_Op& op) const;
        template<typename T>
        std::vector<T> AllReduce(const std::vector<T>& vals, const MPI_Op& op) const;

        template<typename T>
        T Reduce(const T& val, const MPI_Op& op, const int root) const;
        template<typename T>
        std::vector<T> Reduce(const std::vector<T>& vals, const MPI_Op& op, const int root) const;

        template<typename T>
        std::vector<T> Gather(const T& val, const int root) const;

        template<typename T>
        std::vector<T> AllGather(const T& val) const;

        template<typename T>
        std::vector<T> AllToAll(const std::vector<T>& vals) const;

        template<typename T>
        void Send(const T& val, int dest, int tag = 0) const;
        template<typename T>
        void Send(const std::vector<T>& val, int dest, int tag = 0) const;

        template<typename T>
        void Receive(T& val, int src, int tag = 0, MPI_Status* stat = MPI_STATUS_IGNORE) const;
        template<typename T>
        void Receive(std::vector<T>& val, int src, int tag = 0,
                     MPI_Status* stat = MPI_STATUS_IGNORE) const;

      protected:
        /**
         * Constructor to get data needed from an MPI communicator
         * @param communicator
         */
        MpiCommunicator(MPI_Comm communicator, bool willOwn);

        std::shared_ptr<MPI_Comm> commPtr;
    };

    bool operator==(const MpiCommunicator& comm1, const MpiCommunicator& comm2);
    bool operator!=(const MpiCommunicator& comm1, const MpiCommunicator& comm2);

  }
}

#include "net/MpiCommunicator.hpp"

#endif /* HEMELB_NET_MPICOMMUNICATOR_H */
