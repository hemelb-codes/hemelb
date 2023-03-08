// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/decomposition/BasicDecomposition.h"
#include "geometry/GmyReadResult.h"
#include "geometry/LookupTree.h"
#include "net/mpi.h"

namespace hemelb::geometry::decomposition
{

    BasicDecomposition::BasicDecomposition(const GmyReadResult& geometry,
                                           const net::MpiCommunicator& communicator) :
            geometry(geometry), communicator(communicator)
    {
    }

    std::vector<int> BasicDecomposition::Decompose(octree::LookupTree const& tree,
                                                   std::vector<proc_t>& procAssignedToEachBlock)
    {
        // Root node of tree holds total fluid sites
        auto total_sites = tree.levels[0].sites_per_node[0];
        auto target_sites_per_rank = (total_sites - 1) / communicator.Size() + 1;
        // Those blocks which contain at least one fluid site are here
        auto& nonsolid_block_fluid_site_counts = tree.levels[tree.n_levels].sites_per_node;
        auto const n_nonsolid = nonsolid_block_fluid_site_counts.size();
        std::vector<int> rank_for_block(n_nonsolid);

        // We need to return data organised in GMY file order
        // Initialise output to -1 => SOLID, we will overwrite the non-solid below
        std::fill(procAssignedToEachBlock.begin(), procAssignedToEachBlock.end(), -1);
        // We need to know the octree ID to map to 3D coords.
        auto& block_oct_ids = tree.levels[tree.n_levels].node_ids;
        auto& dims = geometry.GetBlockDimensions();
        auto const strides = BlockLocation{dims[1]*dims[2], dims[2], 1};

        int assigned_rank = 0;
        std::decay_t<decltype(nonsolid_block_fluid_site_counts)>::value_type assigned_sites = 0;
        for (int i = 0; i < n_nonsolid; ++i) {
            if (assigned_rank >= communicator.Size())
                throw Exception() << "Trying to assign block to rank that does not exist";

            // This holds only the fluid blocks in octree order
            rank_for_block[i] = assigned_rank;
            // Need to map back to GMY order
            auto ijk = octree::oct_to_ijk(block_oct_ids[i]);
            auto gmy = Dot(ijk, strides);
            procAssignedToEachBlock[gmy] = assigned_rank;

            assigned_sites += nonsolid_block_fluid_site_counts[i];

            if (assigned_sites >= target_sites_per_rank) {
                // Reached target, assign to next rank now
                assigned_rank += 1;
                assigned_sites = 0;
            }
        }
        return rank_for_block;
    }

    void BasicDecomposition::Validate(std::vector<proc_t>& procAssignedToEachBlock)
    {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

        std::vector<proc_t> procForEachBlockRecv = communicator.AllReduce(procAssignedToEachBlock,
                                                                          MPI_MAX);

        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
            if (procAssignedToEachBlock[block] != procForEachBlockRecv[block])
            {
                log::Logger::Log<log::Critical, log::OnePerCore>("At least one other proc thought block %li should be on proc %li but we locally had it as %li",
                                                                 block,
                                                                 procAssignedToEachBlock[block],
                                                                 procForEachBlockRecv[block]);
            }
        }
    }

}
