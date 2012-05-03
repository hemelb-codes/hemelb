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
#include "geometry/ReadResult.h"
#include "geometry/Decomposition.h"

namespace hemelb
{
  namespace geometry
  {

    class GeometryReader
    {
      public:
        GeometryReader(const bool reserveSteeringCore,
                       const lb::lattices::LatticeInfo&,
                       GeometryReadResult& readResult,
                       reporting::Timers &timings);
        ~GeometryReader();

        void LoadAndDecompose(const std::string& dataFilePath);

      private:
        void ReadPreamble();

        void ReadHeader();

        void BlockDecomposition();

        void DivideBlocks(site_t unassignedBlocks,
                          const proc_t unitCount,
                          std::vector<site_t>& blocksOnEachUnit,
                          std::vector<proc_t>& unitForEachBlock,
                          const std::vector<site_t>& fluidSitesPerBlock);

        void ReadInLocalBlocks(const std::vector<proc_t>& unitForEachBlock, const proc_t localRank);

        void DecideWhichBlocksToRead(std::vector<bool>& readBlock,
                                     const std::vector<proc_t>& unitForEachBlock,
                                     const proc_t localRank);

        /**
         * Reads in a single block and ensures it is distributed to all cores that need it.
         * @param offsetSoFar
         * @param procsWantingThisBlock
         * @param blockNumber
         * @param sites
         * @param compressedBytes
         * @param uncompressedBytes
         * @param neededOnThisRank
         */
        void ReadInBlock(MPI_Offset offsetSoFar,
                         const std::vector<proc_t>& procsWantingThisBlock,
                         const site_t blockNumber,
                         const site_t sites,
                         const unsigned int compressedBytes,
                         const unsigned int uncompressedBytes,
                         const int neededOnThisRank);

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

        void ParseBlock(const site_t block, io::writers::xdr::XdrReader& reader);

        SiteReadResult ParseSite(io::writers::xdr::XdrReader& reader);

        /**
         * Calculates the number of the rank used to read in a given block.
         * Intent is to move this into Decomposition class, which will also handle knowledge of which procs to use for reading, and own the decomposition topology.
         *
         * @param blockNumber
         * @return
         */
        proc_t GetReadingCoreForBlock(site_t blockNumber);

        bool Expand(std::vector<util::Vector3D<site_t> >& edgeBlocks,
                    std::vector<util::Vector3D<site_t> >& expansionBlocks,
                    const std::vector<site_t>& fluidSitesPerBlock,
                    std::vector<bool>& blockAssigned,
                    const proc_t currentUnit,
                    std::vector<proc_t>& unitForEachBlock,
                    site_t &blocksOnCurrentUnit,
                    const site_t blocksPerUnit);

        void OptimiseDomainDecomposition(const std::vector<proc_t>& procForEachBlock);

        void ValidateGraphData(const std::vector<idx_t>& vtxDistribn,
                               idx_t localVertexCount,
                               const std::vector<idx_t>& adjacenciesPerVertex,
                               const std::vector<idx_t>& adjacencies);

        void ValidateAllReadData();

        void ValidateProcForEachBlock();

        site_t GetHeaderLength(site_t blockCount) const;

        void GetSiteDistributionArray(std::vector<idx_t>& vertexDistribn,
                                      const std::vector<proc_t>& procForEachBlock,
                                      const std::vector<site_t>& fluidSitesPerBlock) const;

        void GetFirstSiteIndexOnEachBlock(std::vector<idx_t>& firstSiteIndexPerBlock,
                                          const std::vector<idx_t>& vertexDistribution,
                                          const std::vector<proc_t>& procForEachBlock,
                                          const std::vector<site_t>& fluidSitesPerBlock) const;

        void GetAdjacencyData(std::vector<idx_t>& adjacenciesPerVertex,
                              std::vector<idx_t>& localAdjacencies,
                              const idx_t localVertexCount,
                              const std::vector<proc_t>& procForEachBlock,
                              const std::vector<idx_t>& firstSiteIndexPerBlock) const;

        void CallParmetis(std::vector<idx_t>& partitionVector,
                          idx_t localVertexCount,
                          std::vector<idx_t>& vtxDistribn,
                          std::vector<idx_t>& adjacenciesPerVertex,
                          std::vector<idx_t>& adjacencies);

        idx_t* GetMovesList(std::vector<idx_t>& movesFromEachProc,
                            const std::vector<idx_t>& firstSiteIndexPerBlock,
                            const std::vector<proc_t>& procForEachBlock,
                            const std::vector<site_t>& fluidSitesPerBlock,
                            const std::vector<idx_t>& vtxDistribn,
                            const std::vector<idx_t>& partitionVector);

        void RereadBlocks(const std::vector<idx_t>& movesPerProc, const idx_t* movesList, const std::vector<int>& procForEachBlock);

        void ImplementMoves(const std::vector<proc_t>& procForEachBlock,
                            const std::vector<idx_t>& movesFromEachProc,
                            const idx_t* movesList) const;

        proc_t ConvertTopologyRankToGlobalRank(proc_t topologyRank) const;

        static const proc_t HEADER_READING_RANK = 0;
        static const proc_t READING_GROUP_SIZE = HEMELB_READING_GROUP_SIZE;

        const lb::lattices::LatticeInfo& latticeInfo;
        GeometryReadResult& readingResult;
        MPI_File file;
        MPI_Info fileInfo;
        MPI_Comm topologyComm;
        MPI_Group topologyGroup;
        int topologyRank;
        unsigned int topologySize;
        MPI_Comm currentComm;
        int currentCommRank;
        int currentCommSize;
        bool participateInTopology;

        std::vector<site_t> fluidSitesPerBlock;
        std::vector<unsigned int> bytesPerCompressedBlock;
        std::vector<unsigned int> bytesPerUncompressedBlock;
        std::vector<proc_t> procForEachBlock;

        hemelb::reporting::Timers &timings;
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYREADER_H */
