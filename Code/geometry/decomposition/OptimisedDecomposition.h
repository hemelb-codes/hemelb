// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H
#define HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H

#include <vector>
#include <map>
#include "geometry/GmyReadResult.h"
#include "lb/lattices/LatticeInfo.h"
#include "geometry/ParmetisForward.h"
#include "reporting/Timers.h"
#include "net/MpiCommunicator.h"
#include "geometry/SiteData.h"
#include "geometry/GeometryBlock.h"

namespace hemelb::geometry {
    namespace octree { class LookupTree; }

    namespace decomposition
    {
      // Given an initial basic decomposition done at the block level,
      // with all blocks on a process (plus halo) read into the
      // GmyReadResult, use ParMETIS to optimise this.
      //
      // The result of the optimisation is a description of how to
      // change the sites on this rank: those staying, leaving and
      // arriving (see below for details).
      class OptimisedDecomposition
      {
      public:
          // Constructor actually does the optimisation - collective over comm.
          OptimisedDecomposition(reporting::Timers& timers, net::MpiCommunicator comms,
                                 const GmyReadResult& geometry,
                                 const lb::LatticeInfo& latticeInfo);

          // NOTE! All the sites in staying, leaving and arriving are
          // sorted first by block ID and then by intra-block site ID.

          // Get the sites that aren't moved
          inline SiteVec const& GetStaying() const {
              return staying;
          }

          // Get a map (by destination rank) of sites leaving
          inline MovesMap const& GetLeaving() const {
              return leaving;
          }

          // Get a map (by source rank) of sites arriving
          inline MovesMap const& GetArriving() const {
              return arriving;
          }

      private:
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
           * Compile a list of all the moves that need to be made from this processor.
           */
          MovesMap CompileMoveData();

          reporting::Timers& timers; //! Timers for reporting.
          net::MpiCommunicator comms; //! Communicator
          const GmyReadResult& geometry; //! The geometry being optimised.
          octree::LookupTree const& tree;
          const lb::LatticeInfo& latticeInfo; //! The lattice info to optimise for.
          const std::vector<proc_t>& procForBlockOct; //! The initial MPI process for each block, in OCT layout
          const std::vector<U64>& fluidSitesPerBlockOct; //! The number of fluid sites per block, in OCT layout
          std::vector<idx_t> vtxCountPerProc; //! The number of vertices on each process.
          std::vector<idx_t> vtxDistribn; //! The vertex distribution across participating cores.
          std::vector<idx_t> firstSiteIndexPerBlockOct; //! The global contiguous index of the first fluid site on each block.
          std::map<U64, std::vector<U16>> gmySiteIdForBlockOct; //! A map keyed on block OCT id to an array of all the fluid sites' GMY local site ID.
          std::vector<idx_t> adjacenciesPerVertex; //! The number of adjacencies for each local fluid site
          std::vector<idx_t> vertexWeights; //! The weight of each local fluid site
          std::vector<real_t> vertexCoordinates; //! The coordinates of each local fluid site
          std::vector<idx_t> localAdjacencies; //! The list of adjacent vertex numbers for each local fluid site
          std::vector<idx_t> partitionVector; //! The results of the optimisation -- which core each fluid site should go to.
          // The result of the optimisation: the sites staying,
          // leaving and arriving, the latter two organised by the key
          // being the src/dest rank.
          SiteVec staying;
          MovesMap leaving;
          MovesMap arriving;
      };
    }
}
#endif /* HEMELB_GEOMETRY_DECOMPOSITION_OPTIMISEDDECOMPOSITION_H */
