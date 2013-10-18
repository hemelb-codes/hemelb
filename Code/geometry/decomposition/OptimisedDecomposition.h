// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H

#include <vector>
#include "geometry/Geometry.h"
#include "lb/lattices/LatticeInfo.h"
#include "geometry/ParmetisHeader.h"
#include "reporting/Timers.h"
#include "net/MpiCommunicator.h"
#include "geometry/SiteData.h"
#include "geometry/GeometryBlock.h"

namespace hemelb
{
  namespace geometry
  {
    namespace decomposition
    {
      class OptimisedDecomposition
      {
        public:
          OptimisedDecomposition(reporting::Timers& timers, net::MpiCommunicator& comms,
                                 const Geometry& geometry,
                                 const lb::lattices::LatticeInfo& latticeInfo,
                                 const std::vector<proc_t>& procForEachBlock,
                                 const std::vector<site_t>& fluidSitesPerBlock);

          /**
           * Returns a vector with the number of moves coming from each core
           * @return
           */
          inline const std::vector<idx_t>& GetMovesCountPerCore() const
          {
            return allMoves;
          }

          /**
           * Returns a vector with the list of moves
           * @return
           */
          inline const std::vector<idx_t>& GetMovesList() const
          {
            return movesList;
          }

        private:
          typedef util::Vector3D<site_t> BlockLocation;
          /**
           * Populates the vector of vertex weights with different values for each local site type.
           * This allows ParMETIS to more efficiently decompose the system.
           *
           * @return
           */
          void PopulateVertexWeightData(idx_t localVertexCount);
          /**
           * Populates the vertex distribution array in a ParMetis-compatible way. (off-by-1,
           * cumulative count)
           *
           * @param localVertexCount The number of local vertices
           */
          void PopulateSiteDistribution();

          /**
           * Calculate the array of contiguous indices of the first fluid site on each block
           */
          void PopulateFirstSiteIndexOnEachBlock();

          /**
           * Gets the list of adjacencies and the count of adjacencies per local fluid site
           * in a format suitable for use with ParMetis.
           *
           * @param localVertexCount The number of local vertices
           */
          void PopulateAdjacencyData(idx_t localVertexCount);

          /**
           * Perform the call to ParMetis. Returns the result in the partition vector, other
           * parameters are input only. These can't be made const because of the API to ParMetis
           *
           * @param localVertexCount [in] The number of local fluid sites
           */
          void CallParmetis(idx_t localVertexCount);

          /**
           * Populate the list of moves from each proc that we need locally, using the
           * partition vector.
           */
          void PopulateMovesList();

          /**
           * Return true if we should validate.
           * @return
           */
          bool ShouldValidate() const;

          /**
           * Validates the vertex distribution array.
           */
          void ValidateVertexDistribution();

          /**
           * Validates the firstSiteIndexOnEachBlock array
           */
          void ValidateFirstSiteIndexOnEachBlock();

          /**
           * Validate the adjacency data.
           */
          void ValidateAdjacencyData(idx_t localVertexCount);

          /**
           * Sends the adjacency data to the process of lower rank of the two. THIS IS INEFFICIENT.
           * We only do it for validation purposes.
           *
           * @param neighbouringProc
           * @param neighboursAdjacencyCount
           * @param neighboursAdjacencyData Array to receive neighbour's expectations about adjacencies
           * @param expectedAdjacencyData Adjacency data as this core expects it to be
           */
          void SendAdjacencyDataToLowerRankedProc(
              proc_t neighbouringProc, idx_t& neighboursAdjacencyCount,
              std::vector<idx_t>& neighboursAdjacencyData,
              std::multimap<idx_t, idx_t>& expectedAdjacencyData);

          /**
           * Compares this core's and a neighbouring core's version of the adjacency data between
           * them.
           *
           * @param neighbouringProc
           * @param neighboursAdjacencyCount
           * @param neighboursAdjacencyData Array to receive neighbour's expectations about adjacencies
           * @param expectedAdjacencyData Adjacency data as this core expects it to be
           */
          void CompareAdjacencyData(proc_t neighbouringProc, idx_t neighboursAdjacencyCount,
                                    const std::vector<idx_t>& neighboursAdjacencyData,
                                    std::multimap<idx_t, idx_t>& expectedAdjacencyData);

          /**
           * Compile a list of all the moves that need to be made from this processor.
           * @param blockIdLookupByLastSiteIndex
           * @return
           */
          std::vector<idx_t> CompileMoveData(
              std::map<site_t, site_t>& blockIdLookupByLastSiteIndex);

          /**
           * Force some other cores to take info on blocks they might not know they need to know
           * about.
           * @param moveData
           * @param blockIdsIRequireFromX
           */
          void ForceSomeBlocksOnOtherCores(
              std::vector<idx_t>& moveData,
              std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX);

          /**
           * Get the blocks required from every other processor.
           *
           * @param numberOfBlocksRequiredFrom
           * @param blockIdsIRequireFromX
           * @param numberOfBlocksXRequiresFromMe
           * @param blockIdsXRequiresFromMe
           */
          void GetBlockRequirements(
              std::vector<site_t>& numberOfBlocksRequiredFrom,
              std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX,
              std::vector<site_t>& numberOfBlocksXRequiresFromMe,
              std::map<proc_t, std::vector<site_t> >& blockIdsXRequiresFromMe);

          /**
           * Share the number of moves to be made between each pair of processors that need to move
           * data.
           * @param movesForEachLocalBlock
           * @param blockIdsXRequiresFromMe
           * @param coresInterestedInEachBlock
           * @param moveData
           * @param moveDataForEachBlock
           * @param blockIdsIRequireFromX
           * @param movesForEachBlockWeCareAbout
           */
          void ShareMoveCounts(std::map<site_t, idx_t>& movesForEachLocalBlock,
                               std::map<proc_t, std::vector<site_t> >& blockIdsXRequiresFromMe,
                               std::map<site_t, std::vector<proc_t> >& coresInterestedInEachBlock,
                               std::vector<idx_t>& moveData,
                               std::map<site_t, std::vector<idx_t> >& moveDataForEachBlock,
                               std::map<proc_t, std::vector<site_t> >& blockIdsIRequireFromX,
                               std::vector<idx_t>& movesForEachBlockWeCareAbout);

          /**
           * Share the move data between cores
           *
           * @param movesForEachBlockWeCareAbout
           * @param blockIdsIRequireFromX
           * @param blockIdsXRequiresFromMe
           * @param moveDataForEachBlock
           */
          void ShareMoveData(std::vector<idx_t> movesForEachBlockWeCareAbout,
                             std::map<proc_t, std::vector<site_t> > blockIdsIRequireFromX,
                             std::map<proc_t, std::vector<site_t> > blockIdsXRequiresFromMe,
                             std::map<site_t, std::vector<idx_t> > moveDataForEachBlock);

          reporting::Timers& timers; //! Timers for reporting.
          net::MpiCommunicator& comms; //! Communicator
          const Geometry& geometry; //! The geometry being optimised.
          const lb::lattices::LatticeInfo& latticeInfo; //! The lattice info to optimise for.
          const std::vector<proc_t>& procForEachBlock; //! The processor assigned to each block at the moment
          const std::vector<site_t>& fluidSitesPerBlock; //! The number of fluid sites on each block.
          std::vector<idx_t> vtxDistribn; //! The vertex distribution across participating cores.
          std::vector<idx_t> firstSiteIndexPerBlock; //! The global contiguous index of the first fluid site on each block.
          std::vector<idx_t> adjacenciesPerVertex; //! The number of adjacencies for each local fluid site
          std::vector<idx_t> vertexWeights; //! The weight of each local fluid site
          std::vector<idx_t> localAdjacencies; //! The list of adjacent vertex numbers for each local fluid site
          std::vector<idx_t> partitionVector; //! The results of the optimisation -- which core each fluid site should go to.
          std::vector<idx_t> allMoves; //! The list of move counts from each core
          std::vector<idx_t> movesList;
      };
    }
  }
}

#endif /* HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H */
