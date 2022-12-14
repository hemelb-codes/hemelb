// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_DECOMPOSITION_BASICDECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_BASICDECOMPOSITION_H

#include "geometry/GmyReadResult.h"
#include "lb/lattices/LatticeInfo.h"
#include "net/MpiCommunicator.h"
#include "units.h"
#include "util/Vector3D.h"

namespace hemelb::geometry::decomposition
{
      class BasicDecomposition
      {
        public:
          typedef util::Vector3D<site_t> BlockLocation;

          /**
           * Constructor to populate all fields necessary for a decomposition
           *
           * NOTE: We need the geometry and fluidSitesOnEachBlock in order to try to keep
           * contiguous blocks together, and to skip blocks with no fluid sites.
           *
           * @param geometry
           * @param communicator
           * @param fluidSitesOnEachBlock
           */
          BasicDecomposition(const GmyReadResult& geometry, const lb::lattices::LatticeInfo& latticeInfo,
                             const net::MpiCommunicator& communicator,
                             const std::vector<site_t>& fluidSitesOnEachBlock);

          /**
           * Does a basic decomposition of the geometry without requiring any communication;
           * produces a vector of the processor assigned to each block.
           *
           * To make this fast and remove a need to read in the site info about blocks,
           * we assume that neighbouring blocks with any fluid sites on have lattice links
           * between them.
           *
           *  NOTE that the old version of this code used to cope with running on multiple machines,
           * by decomposing fluid sites over machines, then decomposing over processors within one machine.
           * To achieve this here, one could use parmetis's "nparts" parameters to first decompose over machines, then
           * to decompose within machines.
           *
           * @param rankForEachBlock A vector with the processor rank each block has been assigned to.
           */
          void Decompose(std::vector<proc_t>& procAssignedToEachBlock);

          /**
           * Validates that all cores have the same beliefs about which proc is to be assigned
           * to each proc by this decomposition.
           *
           * @param procAssignedToEachBlock This core's decomposition result.
           */
          void Validate(std::vector<proc_t>& procAssignedToEachBlock);

        private:

          /**
           * Does the work of dividing blocks up between processors.
           *
           * The algorithm iterates over processes.
           * We start by assigning the next unassigned block to the current process, then growing out
           * the region by adding new blocks, until the current process has approximately the right
           * number of blocks for the given numbers of blocks / processes. When adding blocks, we prefer
           * blocks that are neighbours of blocks already assigned to the current process.
           *
           * @param processForEachBlock [out] The processor id for each block
           * @param unassignedBlocks [in] The number of blocks yet to be assigned a processor
           * @param geometry [in] The geometry we're decomposing
           * @param processCount [in] The total number of processors
           * @param fluidSitesPerBlock [in] The number of fluid sites in each block
           */
          void DivideBlocks(std::vector<proc_t>& processForEachBlock, site_t unassignedBlocks,
                            const GmyReadResult& geometry, const proc_t processCount,
                            const std::vector<site_t>& fluidSitesPerBlock);

          /**
           * Attempt to expand an already connected volume of blocks assigned to one processor
           * to add additional blocks, already connected to the volume.
           *
           * Returns true if the region was expanded.
           *
           * @param edgeBlocks
           * @param expansionBlocks
           * @param blockAssigned
           * @param currentProcess
           * @param processForEachBlock
           * @param blocksPerProcess
           * @return Returns true if the region was expanded.
           */
          bool Expand(std::vector<BlockLocation>& expansionBlocks, std::vector<bool>& blockAssigned,
                      std::vector<proc_t>& processForEachBlock, site_t &blocksOnCurrentProcess,
                      const std::vector<BlockLocation>& edgeBlocks, const proc_t currentProcess,
                      const site_t blocksPerProcess);

          const GmyReadResult& geometry; //! The geometry being decomposed.
          const lb::lattices::LatticeInfo& latticeInfo; //! The lattice to decompose for.
          const net::MpiCommunicator& communicator; //! The communicator object being decomposed over.
          const std::vector<site_t>& fluidSitesOnEachBlock; //! The number of fluid sites on each block in the geometry.
      };
}
#endif
