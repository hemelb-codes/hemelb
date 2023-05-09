// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "redblood/parallel/GraphBasedCommunication.h"

#include "net/MpiCommunicator.h"
#include "reporting/Timers.h"
#include "util/Iterator.h"

namespace hemelb::redblood::parallel
{

    GlobalCoordsToProcMap ComputeGlobalCoordsToProcMap(net::MpiCommunicator const &comm,
                                                       const geometry::Domain& domain)
    {
        // If we are going to use the octree ordering as-implemented, then
        // sites must have coordinates that fit in 16 bits (per dimension).
        // This can be relaxed of course but requires work. Here we check that.
        {
            auto &largest_coord = domain.GetGlobalSiteMaxes();
            if (std::any_of(largest_coord.begin(), largest_coord.end(),
                            [](LatticeCoordinate x) {
                                return x > LatticeCoordinate(std::numeric_limits<U16>::max());
                            }
            ))
                throw (Exception() << "The domain dimensions in sites do not fit in 16 bits. "
                                      "See comments in " __FILE__ " near line " << __LINE__);
        }

        GlobalCoordsToProcMap coordsToProcMap;
        // Populate map with coordinates of locally owned lattice sites first
        std::vector<Vec16> locallyOwnedSites;
        locallyOwnedSites.reserve(domain.GetLocalFluidSiteCount());
        auto myRank = comm.Rank();
        for (site_t localSiteId = 0; localSiteId < domain.GetLocalFluidSiteCount(); ++localSiteId)
        {
            auto const& globalSiteCoords = domain.GetSite(localSiteId).GetGlobalSiteCoords().as<U16>();
            locallyOwnedSites.push_back(globalSiteCoords);
            coordsToProcMap[globalSiteCoords] = myRank;
        }

        // Exchange coordinates of locally owned lattice sites with neighbours in comms graph
        auto const neighbouringProcs = comm.GetNeighbors();
        if (!neighbouringProcs.empty())
        {
            auto neighSites = comm.AllNeighGatherV(locallyOwnedSites);
            if (std::ssize(neighSites) != comm.GetNeighborsCount())
                throw (Exception() << "Something wrong with neighbourhood");

            // Finish populating map with knowledge coming from neighbours
            for (auto&& [i, p]: util::enumerate(neighbouringProcs)) {
                for (auto const& globalCoord: neighSites[i]) {
                    // lattice sites are uniquely owned, so no chance of coordinates being repeated across processes
                    assert(coordsToProcMap.count(globalCoord) == 0);
                    coordsToProcMap[globalCoord] = p;
                }
            }
        }

        return coordsToProcMap;
    }

    //! \brief Compute ranks this one may communicate with based on which sites are with the RBCs effective size of an edge site.
    std::vector<int> ComputeProcessorNeighbourhood(net::MpiCommunicator const &comm,
                                                   geometry::Domain &domain,
                                                   LatticeDistance cellsEffectiveSize)
    {
        // Keep this sorted for easy lookup
        std::vector<int> ans;

        // Neighbour site IDs that we have checked - kept sorted
        std::vector<U64> checked_ids;
        auto const grid_size = int(std::ceil(cellsEffectiveSize));

        auto valid = [&domain](LatticeVector const& v) {
            return v.IsInRange(domain.GetGlobalSiteMins(), domain.GetGlobalSiteMaxes());
        };
        auto rank = comm.Rank();

        // Loop over this process's edge
        auto first_edge_site_i = domain.GetMidDomainSiteCount();
        auto last_edge_site_i = first_edge_site_i + domain.GetDomainEdgeSiteCount();
        for (auto siteIndex = first_edge_site_i;
             siteIndex < last_edge_site_i;
             ++siteIndex) {
            auto edge_site = domain.GetSite(siteIndex);
            for (int dx = -grid_size; dx <= grid_size; ++dx)
                for (int dy = -grid_size; dy <= grid_size; ++dy)
                    for (int dz = -grid_size; dz <= grid_size; ++dz) {
                        auto delta = util::Vector3D<int>{dx, dy, dz};
                        // Points outside the sphere
                        if (delta.GetMagnitudeSquared() > cellsEffectiveSize*cellsEffectiveSize)
                            continue;
                        auto& site = edge_site.GetGlobalSiteCoords();

                        // Sites with a negative or >= box max index will mess up ID calculation
                        auto neigh = site + delta;
                        if (!valid(neigh))
                            continue;

                        // Have we checked this site before? Search the sorted vector to see.
                        auto neigh_idx = geometry::octree::ijk_to_oct(neigh.as<U16>());
                        //auto neigh_idx = coord_to_id(neigh);
                        auto iter = std::lower_bound(checked_ids.begin(), checked_ids.end(), neigh_idx);
                        if (iter == checked_ids.end()) {
                            // No elem greater than or equal.
                            // Therefore we've not seen it AND can just push onto the end.
                            checked_ids.push_back(neigh_idx);
                        } else if (*iter == neigh_idx) {
                            // We have seen it, nothing else to do.
                            continue;
                        } else {
                            // Not seen it, but have the iterator to the first greater. Insert here.
                            checked_ids.insert(iter, neigh_idx);
                        }

                        // Since we didn't continue, we have an unchecked coord that could have a fluid site.
                        auto [neigh_rank, neigh_local_idx] = domain.GetRankIndexFromGlobalCoords(neigh);
                        if (neigh_rank == SITE_OR_BLOCK_SOLID)
                            // It was not a fluid site
                            continue;

                        if (neigh_rank != rank) {
                            // And it doesn't live on this process.
                            // Find out if the site is edge-of-domain
                            if (domain.IsSiteDomainEdge(neigh_rank, neigh_local_idx)) {
                                // Add that rank to the answer if not already there.
                                auto p_iter = std::lower_bound(ans.begin(), ans.end(), neigh_rank);

                                if (p_iter == ans.end()) {
                                    ans.push_back(neigh_rank);
                                } else if (*p_iter != neigh_rank) {
                                    ans.insert(p_iter, neigh_rank);
                                }
                            }
                        }
                    }
        }
        return ans;
    }

    LatticeDistance ComputeCellsEffectiveSize(TemplateCellContainer const& cellTemplates)
    {
        double maxCellRadius = std::numeric_limits<LatticeDistance>::min();

        for (auto& cellTemplate : cellTemplates)
        {
            maxCellRadius = std::max(maxCellRadius, cellTemplate.second->GetScale());
        }

        return MAXIMUM_SIZE_TO_RADIUS_RATIO * maxCellRadius;
    }

    bool check_neighbourhood_consistency(net::MpiCommunicator const& comm, std::vector<int> const& this_ranks_neighs) {
        // Send to all the ranks I think I have as neighbours
        int const N = std::ssize(this_ranks_neighs);
        std::vector<MPI_Request> sends(N);
        int const rank = comm.Rank();
        int const tag = 4671;
        for (int i = 0; i < N; ++i) {
            HEMELB_MPI_CALL(MPI_Issend,
                            (&rank, 1, MPI_INT, this_ranks_neighs[i], tag, comm, &sends[i]));
        }

        // Signal that I have finished my sends
        MPI_Request barrier = MPI_REQUEST_NULL;
        int my_sends_received = 0;

        // Put the received messages here
        std::vector<int> neighs_sent_to_me;

        int all_sends_received = 0;
        while (!all_sends_received) {
            // Check for arriving messages
            MPI_Status status;
            int msg_available;
            HEMELB_MPI_CALL(MPI_Iprobe,
                            (MPI_ANY_SOURCE, tag, comm, &msg_available, &status));
            if (msg_available) {
                // We have one - check it's correct and push data onto vec
                int nrecv;
                HEMELB_MPI_CALL(MPI_Get_count,
                                (&status, MPI_INT, &nrecv));
                if (nrecv != 1)
                    throw (Exception() << "wrong size");
                int data;
                HEMELB_MPI_CALL(MPI_Recv,
                                (&data, 1, MPI_INT, status.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE));
                if (data != status.MPI_SOURCE)
                    throw (Exception() << "wrong data received");
                neighs_sent_to_me.push_back(data);
            }

            if (!my_sends_received) {
                // check to see if all this rank's sends have been received
                HEMELB_MPI_CALL(MPI_Testall,
                                (N, sends.data(), &my_sends_received, MPI_STATUSES_IGNORE));
                if (my_sends_received)
                    HEMELB_MPI_CALL(MPI_Ibarrier,
                                    (comm, &barrier));
            } else {
                HEMELB_MPI_CALL(MPI_Test,
                                (&barrier, &all_sends_received, MPI_STATUS_IGNORE));
            }
        }

        std::sort(neighs_sent_to_me.begin(), neighs_sent_to_me.end());
        if (this_ranks_neighs.size() != neighs_sent_to_me.size())
            throw (Exception() << "sent and received sizes don't match");

        for (int i = 0; i < N; ++i) {
            if (this_ranks_neighs[i] != neighs_sent_to_me[i])
                throw (Exception() << "Differing neighbour lists");
        }
        return true;
    }

    //! \brief Generates a graph communicator describing the data dependencies for interpolation and spreading
    net::MpiCommunicator CreateGraphComm(net::MpiCommunicator const &comm,
                                         geometry::Domain &domain,
                                         std::shared_ptr<TemplateCellContainer> cellTemplates,
                                         hemelb::reporting::Timers &timings)
    {
        timings[hemelb::reporting::Timers::graphComm].Start();
        auto ranks_i_may_communicate_with = ComputeProcessorNeighbourhood(comm,
                                                                          domain,
                                                                          ComputeCellsEffectiveSize(*cellTemplates));
        check_neighbourhood_consistency(comm, ranks_i_may_communicate_with);
        auto graphComm = comm.DistGraphAdjacent(ranks_i_may_communicate_with);
        timings[hemelb::reporting::Timers::graphComm].Stop();

        return graphComm;
    }
}
