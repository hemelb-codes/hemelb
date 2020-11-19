// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPICOMMUNICATOR_H
#define HEMELB_NET_MPICOMMUNICATOR_H

#include <vector>
#include <map>
#include "net/MpiError.h"
#include <memory>
#include <cassert>

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
         * Copy Constructor
         */
        MpiCommunicator(MpiCommunicator const & comm);

        /**
         * Move Constructor
         */
        MpiCommunicator(MpiCommunicator && comm);

        /**
         * Class has virtual methods so should have virtual d'tor.
         */
        virtual ~MpiCommunicator();

        /**
         * Returns the local rank on the communicator
         * @return
         */
        inline int Rank() const
        {
          assert(localRankInCommunicator >= 0);
          return localRankInCommunicator;
        }

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        inline int Size() const
        {
          assert(communicatorSize > 0);
          return communicatorSize;
        }

        /**
         * Creates a new communicator - see MPI_COMM_GROUP
         * @param Group which is a subset of the group of this communicator.
         * @return New communicator.
         */
        MpiCommunicator Create(const MpiGroup& grp) const;

        //! Copy assignment
        void operator=(MpiCommunicator const &comm);
        //! Move assignment
        void operator=(MpiCommunicator &&comm);

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

        //! \brief Specialization for a vector of variable size
        //! \note Two collective MPI operations are made here, first to get the sizes, then to get
        //! the values.
        template<typename T>
        std::vector<T> Gather(const std::vector<T>& val, const int root) const;

        template <typename T>
        T Scatter(const std::vector<T>& vals, const int root) const;
        template <typename T>
        std::vector<T> Scatter(const std::vector<T>& vals, const size_t n, const int root) const;

        template <typename T>
        std::vector<T> AllGather(const T& val) const;

        /**
         * Performs an all gather operation of fixed size among the neighbours defined in a MPI graph communicator
         * @param val local contribution to all gather operation
         * @return vector with contributions from each neighbour. Use GetNeighbors() to map zero-based indices of the vector to MPI ranks
         */
        template<typename T>
        std::vector<T> AllNeighGather(const T& val) const;

        /**
         * Performs an all gather operation with vectors of variable size among the neighbours defined in a MPI graph communicator
         * @param val vector with local contribution to all gather operation
         * @return vector of vectors with contributions from each neighbour. Use GetNeighbors() to map zero-based indices of outermost vector to MPI ranks
         */
        template<typename T>
        std::vector<std::vector<T>> AllNeighGatherV(const std::vector<T>& val) const;

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

        //! \brief Creates a graph communicator
        //! \param edges: [ [vertices connected to 0], [vertices connected to 1], ...]
        //! \param reorder: Whether nodes can be re-ordered
        MpiCommunicator Graph(std::vector<std::vector<int>> edges, bool reorder = true) const;

        //! \brief Returns graph neighborhood for calling process
        //! \details This communicator must have been created with graph
        std::vector<int> GetNeighbors() const;
        //! \brief Number of neighbors for calling process in a graph communicator
        int GetNeighborsCount() const;

        //! \brief Returns graph neighborhood for a given rank
        //! \details This communicator must have been created with graph
        //! \param whoseNeighbors: process rank whose neighbors we want to know about
        std::vector<int> GetNeighbors(int whoseNeighbors) const;
        //! \brief Number of neighbors in a graph communicator for a given rank
        //! \param whoseNeighbors: process rank whose number of neighbors we want to know about
        int GetNeighborsCount(int whoseNeighbors) const;

        //! A map from the ranks of this communicator to another
        std::map<int, int> RankMap(MpiCommunicator const &valueComm) const;

        //! Splits a communicator
        MpiCommunicator Split(int color, int rank) const;
        //! Splits a communicator
        MpiCommunicator Split(int color) const
        {
          return Split(color, Rank());
        }

        void Barrier() const;

      protected:
        /**
         * Constructor to get data needed from an MPI communicator
         * @param communicator
         */
        MpiCommunicator(MPI_Comm communicator, bool willOwn);

        /**
         * Constructor used during testing to mock a communicator configuration
         * @param localRankInCommunicator hypothetical local rank in communicator
         * @param communicatorSize hypothetical number of MPI processes in communicator
         */
        MpiCommunicator(int localRankInCommunicator, int communicatorSize);

        std::shared_ptr<MPI_Comm> commPtr;

      private:
        //! Size of underlying communicator. Cached for performance
        int communicatorSize;
        //! Local MPI rank on the underlying communicator. Cached for performance
        int localRankInCommunicator;
    };

    bool operator==(const MpiCommunicator& comm1, const MpiCommunicator& comm2);
    bool operator!=(const MpiCommunicator& comm1, const MpiCommunicator& comm2);
  }
}

#include "net/MpiCommunicator.hpp"

#endif /* HEMELB_NET_MPICOMMUNICATOR_H */
