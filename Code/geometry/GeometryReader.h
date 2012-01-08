#ifndef HEMELB_GEOMETRY_GEOMETRYREADER_H
#define HEMELB_GEOMETRY_GEOMETRYREADER_H

#include <vector>
#include <string>

#include "D3Q15.h"
#include "geometry/SiteData.h"
#include "io/writers/xdr/XdrReader.h"
#include "lb/LbmParameters.h"
#include "mpiInclude.h"
#include "parmetis.h"
#include "reporting/Timers.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace geometry
  {
    struct SiteReadResult
    {
      public:
        SiteReadResult() :
            targetProcessor(-1), siteData(BIG_NUMBER3), ioletNormal(-1.), ioletDistance(-1.), wallNormal(NO_VALUE), wallDistance(-1.)
        {
          for (Direction direction = 0; direction < D3Q15::NUMVECTORS - 1; direction++)
          {
            cutDistance[direction] = -1.0;
          }
        }

        // Processor on which to perform lattice-Boltzmann for the site.
        proc_t targetProcessor;

        // Sitedata.
        SiteData siteData;

        // Estimated normal and distance to an inlet or outlet.
        util::Vector3D<double> ioletNormal;
        double ioletDistance;

        // Estimated wall normal (if the site is close to the wall)
        // and distance to wall.
        util::Vector3D<double> wallNormal;
        double wallDistance;

        // cut distances along the non-zero lattice vectors;
        // each one is between 0 and 1 if the surface cuts the corresponding
        // vector or is equal to "NO_VALUE" otherwise
        double cutDistance[D3Q15::NUMVECTORS - 1];
    };

    struct BlockReadResult
    {
      public:
        std::vector<SiteReadResult> Sites;
    };

    struct GeometryReadResult
    {
      public:
        unsigned int stressType;
        util::Vector3D<site_t> blocks;
        site_t blockSize;
        double voxelSize;
        util::Vector3D<double> origin;

        std::vector<BlockReadResult> Blocks;

        site_t GetBlockCount() const
        {
          return blocks.x * blocks.y * blocks.z;
        }

        site_t GetSitesPerBlock() const
        {
          return util::NumericalFunctions::IntegerPower(blockSize, 3);
        }

        /*
         * This needs to commonised with comparable code in the GlobalLatticeData object.
         */
        site_t GetBlockIdFromBlockCoordinates(site_t blockI, site_t blockJ, site_t blockK) const
        {
          return (blockI * blocks.y + blockJ) * blocks.z + blockK;
        }

        /*
         * TODO: This needs to commonised with comparable code in the GlobalLatticeData object.
         */
        site_t GetSiteIdFromSiteCoordinates(site_t siteI, site_t siteJ, site_t siteK) const
        {
          return (siteI * blockSize + siteJ) * blockSize + siteK;
        }
    };

    class GeometryReader
    {
      public:
        GeometryReader(const bool reserveSteeringCore, GeometryReadResult& readResult);
        ~GeometryReader();

        void LoadAndDecompose(std::string& dataFilePath,
                              lb::LbmParameters* bLbmParams,
                              reporting::Timers &timings);

      private:
        struct BlockLocation
        {
            site_t i, j, k;
        };

        void ReadPreamble(lb::LbmParameters* bParams);

        void ReadHeader(site_t* sitesInEachBlock, unsigned int* bytesUsedByBlockInDataFile);

        void BlockDecomposition(const site_t* fluidSitesPerBlock, proc_t* initialProcForEachBlock);

        void DivideBlocks(site_t unassignedBlocks,
                          const proc_t unitCount,
                          site_t* blocksOnEachUnit,
                          proc_t* unitForEachBlock,
                          const site_t* fluidSitesPerBlock);

        void ReadInLocalBlocks(const site_t* sitesPerBlock,
                               const unsigned int* bytesPerBlock,
                               const proc_t* unitForEachBlock,
                               const proc_t localRank);

        void DecideWhichBlocksToRead(bool* readBlock,
                                     const proc_t* unitForEachBlock,
                                     const proc_t localRank);

        /**
         * Reads in a single block and ensures it is distributed to all cores that need it.
         * @param iGlobLatDat
         * @param offsetSoFar
         * @param buffer
         * @param procsWantingThisBlockBuffer
         * @param blockNumber
         * @param sites
         * @param bytes
         * @param neededOnThisRank
         */
        void ReadInBlock(MPI_Offset offsetSoFar,
                         char* buffer,
                         int* procsWantingThisBlockBuffer,
                         const site_t blockNumber,
                         const site_t sites,
                         const unsigned int bytes,
                         const int neededOnThisRank);

        void ParseBlock(const site_t block, io::writers::xdr::XdrReader& reader);

        SiteReadResult ParseSite(io::writers::xdr::XdrReader& reader);

        /**
         * Calculates the number of the rank used to read in a given block.
         *
         * @param blockNumber
         * @return
         */
        proc_t GetReadingCoreForBlock(site_t blockNumber);

        bool Expand(std::vector<BlockLocation>* edgeBlocks,
                    std::vector<BlockLocation>* expansionBlocks,
                    const site_t* fluidSitesPerBlock,
                    bool* blockAssigned,
                    const proc_t currentUnit,
                    proc_t* unitForEachBlock,
                    site_t &blocksOnCurrentUnit,
                    const site_t blocksPerUnit);

        void OptimiseDomainDecomposition(const site_t* sitesPerBlock,
                                         const unsigned int* bytesPerBlock,
                                         const proc_t* procForEachBlock);

        void ValidateGraphData(idx_t* vtxDistribn,
                               idx_t localVertexCount,
                               idx_t* adjacenciesPerVertex,
                               idx_t* adjacencies);

        void ValidateAllReadData();

        void ValidateProcForEachBlock(proc_t* procForEachBlock);

        site_t GetHeaderLength(site_t blockCount) const;

        void GetSiteDistributionArray(idx_t* vertexDistribn,
                                      const proc_t* procForEachBlock,
                                      const site_t* sitesPerBlock) const;

        void GetFirstSiteIndexOnEachBlock(idx_t* firstSiteIndexPerBlock,
                                          const idx_t* vertexDistribution,
                                          const proc_t* procForEachBlock,
                                          const site_t* sitesPerBlock) const;

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
                            const site_t* sitesPerBlock,
                            const idx_t* vtxDistribn,
                            const idx_t* partitionVector);

        void RereadBlocks(const idx_t* movesPerProc,
                          const idx_t* movesList,
                          const site_t* sitesPerBlock,
                          const unsigned int* bytesPerBlock,
                          const int* procForEachBlock);

        void ImplementMoves(const proc_t* procForEachBlock,
                            const idx_t* movesFromEachProc,
                            const idx_t* movesList) const;

        proc_t ConvertTopologyRankToGlobalRank(proc_t topologyRank) const;

        void CreateFileReadType(MPI_Datatype* dataType,
                                const bool* readBlock,
                                const unsigned int* bytesPerBlock) const;

        // The config file starts with:
        // * 1 unsigned int for stress type
        // * 3 unsigned ints for the number of blocks in the x, y, z directions
        // * 1 unsigned int for the block size (number of sites along one edge of a block)
        // * 1 double for the voxel size
        // * 3 doubles for the world-position of site 0
        static const int preambleBytes = 5 * 4 + 4 * 8;
        static const proc_t HEADER_READING_RANK = 0;
        static const proc_t READING_GROUP_SIZE = 5;

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
    };
  }
}

#endif /* HEMELB_GEOMETRY_GEOMETRYREADER_H */
