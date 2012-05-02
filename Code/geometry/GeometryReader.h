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
        GeometryReader(const bool reserveSteeringCore, const lb::lattices::LatticeInfo&, GeometryReadResult& readResult, reporting::Timers &timings);
        ~GeometryReader();

        void LoadAndDecompose(const std::string& dataFilePath);

      private:
        void ReadPreamble();

        void ReadHeader();

        void BlockDecomposition();

        void DivideBlocks(site_t unassignedBlocks,
                          const proc_t unitCount,
                          site_t* blocksOnEachUnit,
                          proc_t* unitForEachBlock,
                          const site_t* fluidSitesPerBlock);

        void ReadInLocalBlocks(const proc_t* unitForEachBlock, const proc_t localRank);

        void DecideWhichBlocksToRead(std::vector<bool>& readBlock, const proc_t* unitForEachBlock, const proc_t localRank);

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
                    const site_t* fluidSitesPerBlock,
                    bool* blockAssigned,
                    const proc_t currentUnit,
                    proc_t* unitForEachBlock,
                    site_t &blocksOnCurrentUnit,
                    const site_t blocksPerUnit);

        void OptimiseDomainDecomposition(const proc_t* procForEachBlock);

        void ValidateGraphData(idx_t* vtxDistribn,
                               idx_t localVertexCount,
                               idx_t* adjacenciesPerVertex,
                               idx_t* adjacencies);

        void ValidateAllReadData();

        void ValidateProcForEachBlock();

        site_t GetHeaderLength(site_t blockCount) const;

        void GetSiteDistributionArray(idx_t* vertexDistribn,
                                      const proc_t* procForEachBlock,
                                      const site_t* fluidSitesPerBlock) const;

        void GetFirstSiteIndexOnEachBlock(idx_t* firstSiteIndexPerBlock,
                                          const idx_t* vertexDistribution,
                                          const proc_t* procForEachBlock,
                                          const site_t* fluidSitesPerBlock) const;

        void GetAdjacencyData(idx_t* adjacenciesPerVertex,
                              idx_t* &localAdjacencies,
                              const idx_t localVertexCount,
                              const proc_t* procForEachBlock,
                              const idx_t* firstSiteIndexPerBlock) const;

        void CallParmetis(idx_t* partitionVector,
                          idx_t localVertexCount,
                          idx_t* vtxDistribn,
                          idx_t* adjacenciesPerVertex,
                          idx_t* adjacencies);

        idx_t* GetMovesList(idx_t* movesFromEachProc,
                            const idx_t* firstSiteIndexPerBlock,
                            const proc_t* procForEachBlock,
                            const site_t* fluidSitesPerBlock,
                            const idx_t* vtxDistribn,
                            const idx_t* partitionVector);

        void RereadBlocks(const idx_t* movesPerProc, const idx_t* movesList, const int* procForEachBlock);

        void ImplementMoves(const proc_t* procForEachBlock,
                            const idx_t* movesFromEachProc,
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
        site_t* fluidSitesPerBlock;
        unsigned int* bytesPerCompressedBlock;
        unsigned int* bytesPerUncompressedBlock;
        proc_t* procForEachBlock;
        hemelb::reporting::Timers &timings;
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYREADER_H */
