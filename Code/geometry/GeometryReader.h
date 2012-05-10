#ifndef HEMELB_GEOMETRY_GEOMETRYREADER_H
#define HEMELB_GEOMETRY_GEOMETRYREADER_H

#include <vector>
#include <string>

#include "io/writers/xdr/XdrReader.h"
#include "lb/lattices/LatticeInfo.h"
#include "lb/LbmParameters.h"
#include "mpiInclude.h"
#include "parmetis.h"
#include "reporting/Timers.h"
#include "util/Vector3D.h"
#include "units.h"
#include "geometry/Geometry.h"
#include "geometry/Decomposition.h"

namespace hemelb
{
  namespace geometry
  {

    class GeometryReader
    {
      public:
        typedef typename util::Vector3D<site_t> BlockLocation;

        GeometryReader(const bool reserveSteeringCore, const lb::lattices::LatticeInfo&, reporting::Timers &timings);
        ~GeometryReader();

        Geometry LoadAndDecompose(const std::string& dataFilePath);

      private:
        Geometry ReadPreamble();

        void ReadHeader(site_t blockCount);

        void BlockDecomposition(const Geometry& geometry);

        /**
         * Get an initial base-level decomposition of the domain macro-blocks over processors (or some other
         * unit of computation, like an entire machine). This will later be improved upon by ParMetis.
         *
         * NOTE: We need the global lattice data and fluid sites per block in order to try to keep
         * contiguous blocks together, and to skip blocks with no fluid sites.
         *
         * @param unassignedBlocks The number of blocks still unassigned
         * @param geometry The geometry object which will end up containing the blocks
         * @param unitCount The number of units (e.g. processors) to divide up the geometry onto
         * @param blocksOnEachUnit The number of blocks assigned to each processor after the decomposition
         * @param unitForEachBlock The rank number for each block after the decomposition
         * @param fluidSitesPerBlock The number of fluid sites on each block
         */
        void DivideBlocks(std::vector<site_t>& blocksOnEachUnit,
                          std::vector<proc_t>& unitForEachBlock,
                          site_t unassignedBlocks,
                          const Geometry& geometry,
                          const proc_t unitCount,
                          const std::vector<site_t>& fluidSitesPerBlock);

        void ReadInBlocksWithHalo(Geometry& geometry,
                                  const std::vector<proc_t>& unitForEachBlock,
                                  const proc_t localRank);

        /**
         * Compile a list of blocks to be read onto this core, including all the ones we perform
         * LB on, and also any of their neighbouring blocks.
         *
         * NOTE: that the skipping-over of blocks without any fluid sites is dealt with by other
         * code.
         *
         * @param geometry [in] Geometry object as it has been read so far
         * @param unitForEachBlock [in] The initial processor assigned to each block
         * @param localRank [in] Local rank number
         * @return Vector with true for each block we should read in.
         */
        std::vector<bool> DecideWhichBlocksToReadIncludingHalo(const Geometry& geometry,
                                                               const std::vector<proc_t>& unitForEachBlock,
                                                               const proc_t localRank);

        /**
         * Reads in a single block and ensures it is distributed to all cores that need it.
         *
         * @param offsetSoFar [in] The offset into the file to read from to get the block.
         * @param geometry [out] The geometry object to populate with info about the block.
         * @param procsWantingThisBlock [in] A list of proc ids where info about this block is required.
         * @param blockNumber [in] The id of the block we're reading.
         * @param neededOnThisRank [in] A boolean indicating whether the block is required locally.
         */
        void ReadInBlock(MPI_Offset offsetSoFar,
                         Geometry& geometry,
                         const std::vector<proc_t>& procsWantingThisBlock,
                         const site_t blockNumber,
                         const bool neededOnThisRank);

        /**
         * Decompress the block data. Uses the known number of sites to get an
         * upper bound on the uncompressed data to simplify the code and avoid
         * reallocation.
         * @param compressed
         * @param sites
         * @return
         */
        std::vector<char> DecompressBlockData(const std::vector<char>& compressed,
                                              const unsigned int uncompressedBytes);

        void ParseBlock(Geometry& geometry, const site_t block, io::writers::xdr::XdrReader& reader);

        /**
         * Parse the next site from the XDR reader. Note that we return by copy here.
         * @param reader
         * @return
         */
        GeometrySite ParseSite(io::writers::xdr::XdrReader& reader);

        /**
         * Calculates the number of the rank used to read in a given block.
         * Intent is to move this into Decomposition class, which will also handle knowledge of which procs to use for reading, and own the decomposition topology.
         *
         * @param blockNumber
         * @return
         */
        proc_t GetReadingCoreForBlock(site_t blockNumber);

        bool Expand(std::vector<BlockLocation>& edgeBlocks,
                    std::vector<BlockLocation>& expansionBlocks,
                    const std::vector<site_t>& fluidSitesPerBlock,
                    std::vector<bool>& blockAssigned,
                    const proc_t currentUnit,
                    std::vector<proc_t>& unitForEachBlock,
                    site_t &blocksOnCurrentUnit,
                    const site_t blocksPerUnit,
                    const Geometry& geometry);

        /**
         * Optimise the domain decomposition using ParMetis. We take this approach because ParMetis
         * is more efficient when given an initial decomposition to start with.
         * @param geometry
         * @param procForEachBlock
         */
        void OptimiseDomainDecomposition(Geometry& geometry, const std::vector<proc_t>& procForEachBlock);

        void ValidateGraphData(const std::vector<idx_t>& vtxDistribn,
                               idx_t localVertexCount,
                               const std::vector<idx_t>& adjacenciesPerVertex,
                               const std::vector<idx_t>& adjacencies);

        void ValidateGeometry(const Geometry& geometry);

        void ValidateProcForEachBlock(site_t blockCount);

        /**
         * Get the length of the header section, given the number of blocks.
         *
         * @param blockCount The number of blocks.
         * @return
         */
        site_t GetHeaderLength(site_t blockCount) const;

        std::vector<idx_t> GetSiteDistributionArray(site_t blockCount,
                                                    const std::vector<proc_t>& procForEachBlock,
                                                    const std::vector<site_t>& fluidSitesPerBlock) const;

        void GetFirstSiteIndexOnEachBlock(std::vector<idx_t>& firstSiteIndexPerBlock,
                                          site_t blockCount,
                                          const std::vector<idx_t>& vertexDistribution,
                                          const std::vector<proc_t>& procForEachBlock,
                                          const std::vector<site_t>& fluidSitesPerBlock) const;

        /**
         * Gets the list of adjacencies and the count of adjacencies per local fluid site
         * in a format suitable for use with ParMetis.
         *
         * @param adjacenciesPerVertex [out] The number of adjacencies for each local fluid site
         * @param localAdjacencies [out] The list of adjacent vertex numbers for each local fluid site
         * @param geometry The geometry that's been read in
         * @param localVertexCount The number of local vertices
         * @param procForEachBlock The initial processor for each block
         * @param firstSiteIndexPerBlock The universal contiguous index of the first fluid site
         * on each block.
         */
        void GetAdjacencyData(std::vector<idx_t>& adjacenciesPerVertex,
                              std::vector<idx_t>& localAdjacencies,
                              const Geometry& geometry,
                              const idx_t localVertexCount,
                              const std::vector<proc_t>& procForEachBlock,
                              const std::vector<idx_t>& firstSiteIndexPerBlock) const;

        /**
         * Perform the call to ParMetis. Returns the result in the partition vector, other
         * parameters are input only. These can't be made const because of the API to ParMetis
         *
         * @param partitionVector [out] The result from ParMetis, giving the destination proc of each
         * local fluid site
         * @param localVertexCount [in] The number of local fluid sites
         * @param vtxDistribn [in] The number of fluid sites on each cores (in ParMetis-compatible
         * format)
         * @param adjacenciesPerVertex [in] The number of adjacencies each local fluid site has, in
         * ParMetis-compatible format
         * @param adjacencies [in] The list of adjacent fluid site ids for each local fluid site in
         * order
         */
        void CallParmetis(std::vector<idx_t>& partitionVector,
                          idx_t localVertexCount,
                          std::vector<idx_t>& vtxDistribn,
                          std::vector<idx_t>& adjacenciesPerVertex,
                          std::vector<idx_t>& adjacencies);

        idx_t* GetMovesList(std::vector<idx_t>& movesFromEachProc,
                            const Geometry& geometry,
                            const std::vector<idx_t>& firstSiteIndexPerBlock,
                            const std::vector<proc_t>& procForEachBlock,
                            const std::vector<site_t>& fluidSitesPerBlock,
                            const std::vector<idx_t>& vtxDistribn,
                            const std::vector<idx_t>& partitionVector);

        void RereadBlocks(Geometry& geometry,
                          const std::vector<idx_t>& movesPerProc,
                          const idx_t* movesList,
                          const std::vector<int>& procForEachBlock);

        void ImplementMoves(Geometry& geometry,
                            const std::vector<proc_t>& procForEachBlock,
                            const std::vector<idx_t>& movesFromEachProc,
                            const idx_t* movesList) const;

        proc_t ConvertTopologyRankToGlobalRank(proc_t topologyRank) const;

        /**
         * True if we should validate the geometry.
         * @return
         */
        bool ShouldValidate() const;

        //! The rank which reads in the header information.
        static const proc_t HEADER_READING_RANK = 0;
        //! The number of cores (0-READING_GROUP_SIZE-1) that read files in parallel
        static const proc_t READING_GROUP_SIZE = HEMELB_READING_GROUP_SIZE;

        //! Info about the connectivity of the lattice.
        const lb::lattices::LatticeInfo& latticeInfo;
        //! File accessed to read in the geometry data.
        MPI_File file;
        //! Information about the file, to give cues and hints to MPI.
        MPI_Info fileInfo;
        //! MPI Communicator for all ranks that will need a slice of the geometry.
        MPI_Comm topologyComm;
        //! MPI Group of the ranks that will take a piece of the geometry.
        MPI_Group topologyGroup;
        //! This core's rank within the topology group.
        int topologyRank;
        //! The size of the topology group.
        unsigned int topologySize;
        // TODO: This was never a good plan, better code design will avoid the need for it.
        //! The communicator currently in use.
        MPI_Comm currentComm;
        //! This core's rank within the current communicator
        int currentCommRank;
        //! The size of the current communicator
        int currentCommSize;
        //! True iff this rank is participating in the domain decomposition.
        bool participateInTopology;

        //! The number of fluid sites on each block in the geometry
        std::vector<site_t> fluidSitesOnEachBlock;
        //! The number of bytes each block takes up while still compressed.
        std::vector<unsigned int> bytesPerCompressedBlock;
        //! The number of bytes each block takes up when uncompressed.
        std::vector<unsigned int> bytesPerUncompressedBlock;
        //! The processor assigned to each block.
        std::vector<proc_t> principalProcForEachBlock;

        //! Timings object for recording the time taken for each step of the domain decomposition.
        hemelb::reporting::Timers &timings;
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYREADER_H */
