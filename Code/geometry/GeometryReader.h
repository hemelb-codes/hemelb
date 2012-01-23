#ifndef HEMELB_GEOMETRY_GEOMETRYREADER_H
#define HEMELB_GEOMETRY_GEOMETRYREADER_H

#include <vector>
#include <string>

#include "D3Q15.h"
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
    //! Model of the data read in about a link from one site in one direction.
    struct LinkReadResult
    {
      public:
        //! Enumeration of the different intersections the link might make between the current
        //! site and the next lattice point in this direction: no intersection,
        //! intersection with a vessel wall and intersection with an inlet or outlet.
        enum IntersectionType
        {
          NO_INTERSECTION = 0,
          WALL_INTERSECTION = 1,
          INLET_INTERSECTION = 2,
          OUTLET_INTERSECTION = 3
        } type;

        //! Default constructor. Has no intersection, nonsense values for intersection distance
        //! and iolet id.
        LinkReadResult() :
            type(NO_INTERSECTION), distanceToIntersection(-1.0), ioletId(-1)
        {
        }

        //! The proportion of the lattice vector traversed before an intersection is hit.
        float distanceToIntersection;

        //! The identifier of the inlet or outlet hit along the lattice vector (if one is hit).
        int ioletId;
    };

    /***
     * Model of the data for a site, as contained within a geometry file.
     * this data will be broken up and placed in various arrays in hemelb::Geometry::LatticeData
     */
    struct SiteReadResult
    {
      public:
        //! Basic constructor for solid and fluid sites.
        SiteReadResult(bool siteIsFluid) :
            targetProcessor(siteIsFluid ?
              -1 :
              BIG_NUMBER2), isFluid(siteIsFluid)
        {
        }

        SiteReadResult(const SiteReadResult& other) :
            targetProcessor(other.targetProcessor), isFluid(other.isFluid), links(other.links)
        {

        }

        //! Processor on which to perform lattice-Boltzmann for the site.
        proc_t targetProcessor;

        //! True iff the site is fluid, i.e. it is within the geometry and we will be doing
        //! lattice-Boltzmann with it.
        bool isFluid;

        //! A vector of the link data for each direction in the lattice currently being used
        //! (NOT necessarily the same as the lattice used by the geometry file).
        std::vector<LinkReadResult> links;
    };

    /***
     * Model of the information stored for a block in a geometry file.
     * Just gives the array of sites
     */
    struct BlockReadResult
    {
      public:
        std::vector<SiteReadResult> Sites;
    };

    /***
     * Model of the information in a geometry file
     */
    struct GeometryReadResult
    {
      public:
        util::Vector3D<site_t> blocks; //! The count of blocks in each direction
        site_t blockSize; //! Size of a block, in sites.
        double voxelSize; //! Size of a block, in real-world units.
        util::Vector3D<double> origin;

        std::vector<BlockReadResult> Blocks; //! Array of Block models

        /***
         * Total count of blocks in the site.
         * @return count of blocks in the site.
         */
        site_t GetBlockCount() const
        {
          return blocks.x * blocks.y * blocks.z;
        }

        site_t GetSitesPerBlock() const
        {
          return util::NumericalFunctions::IntegerPower(blockSize, 3);
        }

        /***
         * Get the i.d. of a block, i.e. the one-d coordinate, from the three-d coordinate.
         * @todo Use this to replace LatticeData::GetBlockIdFromBlockCoords
         */
        site_t GetBlockIdFromBlockCoordinates(site_t blockI, site_t blockJ, site_t blockK) const
        {
          return (blockI * blocks.y + blockJ) * blocks.z + blockK;
        }

        /***
         * Get the i.d. of a site, i.e. the one-d coordinate, from the three-d coordinate.
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

        void LoadAndDecompose(std::string& dataFilePath, reporting::Timers &timings);

      private:
        void ReadPreamble();

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

        void DecideWhichBlocksToRead(bool* readBlock, const proc_t* unitForEachBlock, const proc_t localRank);

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

        bool Expand(std::vector<util::Vector3D<site_t> >& edgeBlocks,
                    std::vector<util::Vector3D<site_t> >& expansionBlocks,
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
