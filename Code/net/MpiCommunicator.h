// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_MPICOMMUNICATOR_H
#define HEMELB_NET_MPICOMMUNICATOR_H

#include <map>
#include <memory>
#include <numeric>
#include <span>
#include <vector>

#include "hassert.h"
#include "net/MpiError.h"

namespace hemelb::net
{
    class MpiGroup;

    // This type hold the data that results from operations like
    // Allgatherv. The data is stored contiguously in the `data`
    // vector and the displacements (an array of int, as required
    // by MPI) in `displacements`. This has a size of the number
    // of participating processes, plus one.
    //
    // This type presents an interface similar to an array of
    // arrays: specifically, indexing returns a span over the
    // corresponding process's data.
    template <typename T>
    struct displaced_data {
        // Element i gives the index into data for the data from rank i
        // Element MPI_Comm_size gives the total size of data
        std::vector<int> displacements;
        std::vector<T> data;

        // Default is empty.
        displaced_data() = default;

        // Construct from a vector of sizes, computing the displacements.
        displaced_data(std::vector<int> const& sizes) {
            auto const P = sizes.size();
            displacements.resize(P + 1);
            displacements[0] = 0;
            std::inclusive_scan(
                    sizes.begin(), sizes.end(),
                    displacements.begin() + 1
            );
            data.resize(displacements[P]);
        }

        // Size should be the number of participating processes from the collective.
        std::size_t size() const {
            return displacements.size() - 1;
        }

        // Return a span over the data belonging to the given process.
        std::span<T> operator[](std::size_t i) {
            HASSERT(i >= 0 && i < displacements.size());
            auto start = displacements[i];
            auto count = displacements[i+1] - start;
            return std::span<T>{data.data() + start, unsigned(count)};
        }
    };

    // Holds an MPI communicator and exposes communication functions
    // via members. It will own the underlying MPI_Comm (i.e. it will
    // call MPI_Comm_free) if created from another MpiCommunicator
    // instance. Normal value semantics, so copyable and movable.
    // Note that copying does not duplicate the communicator, it just
    // increments the reference count.
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
         * Returns the local rank on the communicator
         * @return
         */
        inline int Rank() const
        {
          HASSERT(localRankInCommunicator >= 0);
          return localRankInCommunicator;
        }

        /**
         * Returns the size of the communicator (i.e. total number of procs involved).
         * @return
         */
        inline int Size() const
        {
          HASSERT(communicatorSize > 0);
          return communicatorSize;
        }

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

        // Broadcast overloads - all are collective as use MPI_Bcast
        //
        // Scalar
        template<typename T>
        void Broadcast(T& val, int root) const;
        // Array of scalar (N can be dynamic_extent)
        // Note that non-root processes must have the same size span!!
        template<typename T, std::size_t N>
        void Broadcast(std::span<T, N> vals, int root) const;
        // Special case for string that does length then data.
        void Broadcast(std::string& val, const int root) const;

        template<typename T>
        T AllReduce(const T& val, const MPI_Op& op) const;
        template<typename T>
        void AllReduceInPlace(const std::span<T>& vals, const MPI_Op& op) const;
        template<typename T>
        std::vector<T> AllReduce(const std::vector<T>& vals, const MPI_Op& op) const;

	template <typename T>
	T Scan(const T& val, const MPI_Op& op) const;

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

        template <typename T>
        displaced_data<T> AllGatherV(const std::vector<T>& vals) const;

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
        displaced_data<T> AllNeighGatherV(const std::vector<T>& val) const;

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

        //! \brief Create a distributed graph communicator assuming unweighted and bidirectional communication.
        MpiCommunicator DistGraphAdjacent(std::vector<int> my_neighbours, bool reorder = true) const;

        //! \brief Returns graph neighborhood for calling process
        //! \details This communicator must have been created with graph
        std::vector<int> GetNeighbors() const;
        //! \brief Number of neighbors for calling process in a graph communicator
        int GetNeighborsCount() const;

        //! A map from the ranks of this communicator to another
        std::map<int, int> RankMap(MpiCommunicator const &valueComm) const;

        //! Splits a communicator
        MpiCommunicator Split(int color, int rank) const;
        //! Splits a communicator
        MpiCommunicator Split(int color) const
        {
          return Split(color, Rank());
        }

        MpiCommunicator SplitType(int type = MPI_COMM_TYPE_SHARED) const;

        void Barrier() const;

      protected:
        /**
         * Constructor to get data needed from an MPI communicator
         * @param communicator
         */
        MpiCommunicator(MPI_Comm communicator, bool willOwn);

        struct mock_ctor_tag {};
        /**
         * Constructor used during testing to mock a communicator configuration
         * @param localRankInCommunicator hypothetical local rank in communicator
         * @param communicatorSize hypothetical number of MPI processes in communicator
         */
        MpiCommunicator(mock_ctor_tag, int localRankInCommunicator, int communicatorSize);

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

#include "net/MpiCommunicator.hpp"

#endif /* HEMELB_NET_MPICOMMUNICATOR_H */
