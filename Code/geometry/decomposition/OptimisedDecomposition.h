#ifndef HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H

#include <vector>
#include "geometry/Geometry.h"
#include "lb/lattices/LatticeInfo.h"
#include "parmetis.h"
#include "reporting/Timers.h"
#include "topology/Communicator.h"

namespace hemelb
{
  namespace geometry
  {
    namespace decomposition
    {
      class OptimisedDecomposition
      {
        public:
          OptimisedDecomposition(reporting::Timers& timers,
                                 topology::Communicator& comms,
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
          typedef typename util::Vector3D<site_t> BlockLocation;
          /**
           * Populates the vertex distribution array in a ParMetis-compatible way. (off-by-1,
           * cumulative count)
           *
           * @return
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
          void SendAdjacencyDataToLowerRankedProc(proc_t neighbouringProc,
                                                  idx_t& neighboursAdjacencyCount,
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
          void CompareAdjacencyData(proc_t neighbouringProc,
                                    idx_t neighboursAdjacencyCount,
                                    const std::vector<idx_t>& neighboursAdjacencyData,
                                    std::multimap<idx_t, idx_t>& expectedAdjacencyData);

          reporting::Timers& timers; //! Timers for reporting.
          topology::Communicator& comms; //! Communicator
          const Geometry& geometry; //! The geometry being optimised.
          const lb::lattices::LatticeInfo& latticeInfo; //! The lattice info to optimise for.

          const std::vector<proc_t>& procForEachBlock; //! The processor assigned to each block at the moment
          const std::vector<site_t>& fluidSitesPerBlock; //! The number of fluid sites on each block.

          std::vector<idx_t> vtxDistribn; //! The vertex distribution across participating cores.
          std::vector<idx_t> firstSiteIndexPerBlock; //! The global contiguous index of the first fluid site on each block.

          std::vector<idx_t> adjacenciesPerVertex; //! The number of adjacencies for each local fluid site
          std::vector<idx_t> localAdjacencies; //! The list of adjacent vertex numbers for each local fluid site

          std::vector<idx_t> partitionVector; //! The results of the optimisation -- which core each fluid site should go to.

          std::vector<idx_t> allMoves; //! The list of move counts from each core
          std::vector<idx_t> movesList; //! The list of all moves (in order of source core).
      };
    }
  }
}

#endif /* HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H */
