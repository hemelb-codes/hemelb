// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "geometry/decomposition/BasicDecomposition.h"
#include "net/mpi.h"

namespace hemelb
{
  namespace geometry
  {
    namespace decomposition
    {

      BasicDecomposition::BasicDecomposition(const Geometry& geometry,
                                             const lb::lattices::LatticeInfo& latticeInfo,
                                             const net::MpiCommunicator& communicator,
                                             const std::vector<site_t>& fluidSitesOnEachBlock) :
        geometry(geometry), latticeInfo(latticeInfo), communicator(communicator),
            fluidSitesOnEachBlock(fluidSitesOnEachBlock)
      {
      }

      void BasicDecomposition::Decompose(std::vector<proc_t>& procAssignedToEachBlock)
      {
        // Keep a count of the number of non-empty blocks that haven't yet been assigned
        // a processor.
        site_t unvisitedFluidBlockCount = 0;
        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
          if (fluidSitesOnEachBlock[block] != 0)
          {
            ++unvisitedFluidBlockCount;
          }
        }

        // Divide blocks between the processors.
        procAssignedToEachBlock.resize(geometry.GetBlockCount());

        DivideBlocks(procAssignedToEachBlock,
                     unvisitedFluidBlockCount,
                     geometry,
                     communicator.Size(),
                     fluidSitesOnEachBlock);
      }

      void BasicDecomposition::Validate(std::vector<proc_t>& procAssignedToEachBlock)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

        std::vector<proc_t> procForEachBlockRecv = communicator.AllReduce(procAssignedToEachBlock, MPI_MAX);

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

      void BasicDecomposition::DivideBlocks(std::vector<proc_t>& unitForEachBlock,
                                            site_t unassignedBlocks,
                                            const Geometry& geometry,
                                            const proc_t unitCount,
                                            const std::vector<site_t>& fluidSitesPerBlock)
      {
        // Initialise the unit being assigned to, and the approximate number of blocks
        // required on each unit.
        proc_t currentUnit = 0;

        site_t targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (communicator.Size()));

        // Create an array to monitor whether each block has been assigned yet.
        std::vector<bool> blockAssigned(geometry.GetBlockCount(), false);

        // Create lists of the current edge of blocks on the current proc and the edge being expanded into
        std::vector<BlockLocation> currentEdge;
        std::vector<BlockLocation> expandedEdge;

        site_t blockNumber = -1;

        // Domain Decomposition.  Pick a site. Set it to the rank we are
        // looking at. Find its neighbours and put those on the same
        // rank, then find the next-nearest neighbours, etc. until we
        // have a completely joined region, or there are enough fluid
        // sites on the rank.  In the former case, start again at
        // another site. In the latter case, move on to the next rank.
        // Do this until all sites are assigned to a rank. There is a
        // high chance of of all sites on a rank being joined.

        site_t blocksOnCurrentProc = 0;

        // Iterate over all blocks.
        for (site_t blockCoordI = 0; blockCoordI < geometry.GetBlockDimensions().x; blockCoordI++)
        {
          for (site_t blockCoordJ = 0; blockCoordJ < geometry.GetBlockDimensions().y; blockCoordJ++)
          {
            for (site_t blockCoordK = 0; blockCoordK < geometry.GetBlockDimensions().z; blockCoordK++)
            {
              // Block number is the number of the block we're currently on.
              blockNumber++;

              // If the array of proc rank for each site is nullptr, we're on an all-solid block.
              // Alternatively, if this block has already been assigned, move on.
              if (fluidSitesPerBlock[blockNumber] == 0)
              {
                unitForEachBlock[blockNumber] = -1;
                continue;
              }
              else if (blockAssigned[blockNumber])
              {
                continue;
              }

              // Assign this block to the current unit.
              blockAssigned[blockNumber] = true;
              unitForEachBlock[blockNumber] = currentUnit;

              ++blocksOnCurrentProc;

              // Record the location of this initial site.
              currentEdge.clear();
              BlockLocation lNew(blockCoordI, blockCoordJ, blockCoordK);
              currentEdge.push_back(lNew);

              // The subdomain can grow.
              bool isRegionGrowing = true;

              // While the region can grow (i.e. it is not bounded by solids or visited
              // sites), and we need more sites on this particular rank.
              while (blocksOnCurrentProc < targetBlocksPerUnit && isRegionGrowing)
              {
                expandedEdge.clear();

                // Sites added to the edge of the mClusters during the iteration.
                isRegionGrowing = Expand(expandedEdge,
                                         blockAssigned,
                                         unitForEachBlock,
                                         blocksOnCurrentProc,
                                         currentEdge,
                                         currentUnit,
                                         targetBlocksPerUnit);

                // When the new layer of edge sites has been found, swap the buffers for
                // the current and new layers of edge sites.
                currentEdge.swap(expandedEdge);
              }

              // If we have enough sites, we have finished.
              if (blocksOnCurrentProc >= targetBlocksPerUnit)
              {
                ++currentUnit;

                unassignedBlocks -= blocksOnCurrentProc;
                targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (unitCount - currentUnit));

                blocksOnCurrentProc = 0;
              }
              // If not, we have to start growing a different region for the same rank:
              // region expansions could get trapped.

            } // Block co-ord k
          } // Block co-ord j
        } // Block co-ord i
      }

      bool BasicDecomposition::Expand(std::vector<BlockLocation>& expansionBlocks,
                                      std::vector<bool>& blockAssigned,
                                      std::vector<proc_t>& unitForEachBlock,
                                      site_t &blocksOnCurrentUnit,
                                      const std::vector<BlockLocation>& edgeBlocks,
                                      const proc_t currentUnit,
                                      const site_t blocksPerUnit)
      {
        bool regionExpanded = false;

        // For sites on the edge of the domain (sites_a), deal with the neighbours.
        for (unsigned int edgeBlockId = 0; (edgeBlockId < edgeBlocks.size()) && (blocksOnCurrentUnit < blocksPerUnit); edgeBlockId++)
        {
          const BlockLocation& edgeBlockCoords = edgeBlocks[edgeBlockId];

          for (Direction direction = 1; direction < latticeInfo.GetNumVectors() && blocksOnCurrentUnit < blocksPerUnit; direction++)
          {
            // Record neighbour location.
            BlockLocation neighbourCoords = edgeBlockCoords + latticeInfo.GetVector(direction);

            // Move on if neighbour is outside the bounding box.
            if (!geometry.AreBlockCoordinatesValid(neighbourCoords))
            {
              continue;
            }

            // Move on if the neighbour is in a block of solids (in which case
            // the pointer to ProcessorRankForEachBlockSite is nullptr) or it is solid or has already
            // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
            // was initialized in lbmReadConfig in io.cc.

            site_t neighBlockId = geometry.GetBlockIdFromBlockCoordinates(neighbourCoords.x,
                                                                          neighbourCoords.y,
                                                                          neighbourCoords.z);

            // Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
            if (fluidSitesOnEachBlock[neighBlockId] == 0 || blockAssigned[neighBlockId])
            {
              continue;
            }

            // Set the rank for a neighbour and update the fluid site counters.
            blockAssigned[neighBlockId] = true;
            unitForEachBlock[neighBlockId] = currentUnit;
            ++blocksOnCurrentUnit;

            // Neighbour was found, so the region can grow.
            regionExpanded = true;

            // Record the location of the neighbour.
            expansionBlocks.push_back(neighbourCoords);
          }
        }

        return regionExpanded;
      }

    } /* namespace decomposition */
  } /* namespace geometry */
} /* namespace hemelb */
