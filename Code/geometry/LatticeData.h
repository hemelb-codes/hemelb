#ifndef HEMELB_GEOMETRY_LATTICEDATA_H
#define HEMELB_GEOMETRY_LATTICEDATA_H

#include <cstdio>

#include "net/net.h"
#include "parmetis.h"
#include "D3Q15.h"
#include "constants.h"
#include "configuration/SimConfig.h"
#include "mpiInclude.h"
#include "geometry/BlockTraverser.h"
#include "io/writers/xdr/XdrReader.h"
#include "reporting/Timers.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    class LatticeData
    {
      public:
        enum SiteType
        {
          // These must be consistent with the setup tool
          SOLID_TYPE = 0U,
          FLUID_TYPE = 1U,
          INLET_TYPE = 2U,
          OUTLET_TYPE = 3U
        };

        static LatticeData* Load(const bool reserveSteeringCore,
                                 std::string& dataFilePath,
                                 reporting::Timers &timings);

        virtual ~LatticeData();

        void SwapOldAndNew();
        void SendAndReceive(net::Net* net);
        void CopyReceived();

        const double* GetNormalToWall(site_t iSiteIndex) const;

        site_t GetXSiteCount() const;
        site_t GetYSiteCount() const;
        site_t GetZSiteCount() const;

        site_t GetXBlockCount() const;
        site_t GetYBlockCount() const;
        site_t GetZBlockCount() const;

        distribn_t GetVoxelSize() const;
        distribn_t GetXOrigin() const;
        distribn_t GetYOrigin() const;
        distribn_t GetZOrigin() const;

        unsigned int GetLog2BlockSize() const;

        site_t GetBlockSize() const;
        site_t GetBlockCount() const;

        site_t GetSitesPerBlockVolumeUnit() const;

        site_t GetBlockIdFromBlockCoords(site_t i, site_t j, site_t k) const;

        bool IsValidBlock(site_t i, site_t j, site_t k) const;
        bool IsValidLatticeSite(site_t i, site_t j, site_t k) const;

        const proc_t
        * GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const;

        BlockData* GetBlock(site_t blockNumber) const;
        BlockTraverser GetBlockTraverser() const;

        distribn_t* GetFOld(site_t siteNumber) const;
        distribn_t* GetFNew(site_t siteNumber) const;
        site_t GetLocalFluidSiteCount() const;
        SiteType GetSiteType(site_t iSiteIndex) const;
        int GetBoundaryId(site_t iSiteIndex) const;
        site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const;
        bool HasBoundary(const site_t iSiteIndex, const int iDirection) const;
        double GetCutDistance(site_t iSiteIndex, int iDirection) const;
        unsigned int GetSiteData(site_t iSiteIndex) const;
        unsigned int GetContiguousSiteId(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;
        const util::Vector3D<site_t>
        GetGlobalCoords(site_t blockNumber, const util::Vector3D<site_t>& localSiteCoords) const;

        site_t GetInnerSiteCount() const;
        site_t GetInnerCollisionCount(unsigned int collisionType) const;
        site_t GetInterCollisionCount(unsigned int collisionType) const;
        unsigned int GetCollisionType(unsigned int site_data) const;
        const site_t* GetFluidSiteCountsOnEachProc() const;
        site_t GetFluidSiteCountOnProc(proc_t proc) const;
        site_t GetTotalFluidSites() const;
        const util::Vector3D<site_t>& GetGlobalSiteMins() const;
        const util::Vector3D<site_t>& GetGlobalSiteMaxes() const;

      protected:
        class LocalLatticeData
        {
          public:
            LocalLatticeData();
            LocalLatticeData(site_t iLocalFluidSites);
            ~LocalLatticeData();

            site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const;
            double GetCutDistance(site_t iSiteIndex, int iDirection) const;
            bool HasBoundary(const site_t iSiteIndex, const int iDirection) const;
            int GetBoundaryId(site_t iSiteIndex) const;
            const double *GetNormalToWall(site_t iSiteIndex) const;
            SiteType GetSiteType(site_t iSiteIndex) const;
            site_t GetLocalFluidSiteCount() const;

            void SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue);
            void SetWallNormal(site_t iSiteIndex, const double iNormal[3]);
            void SetDistanceToWall(site_t iSiteIndex,
                                   const double iCutDistance[D3Q15::NUMVECTORS - 1]);

            void SetSharedSiteCount(site_t iSharedCount);

          public:
            site_t my_inner_sites;
            site_t my_inner_collisions[COLLISION_TYPES];
            site_t my_inter_collisions[COLLISION_TYPES];

            distribn_t *FOld;
            distribn_t *FNew;

            // TODO sadly this has to be public, due to some budgetry in the way we determine site type.
            // SiteType || FluidSite and SiteType && FluidSite have different significances...
            unsigned int *mSiteData;

          private:
            site_t LocalFluidSites;
            site_t* mFNeighbours;
            double *mDistanceToWall;
            double *mWallNormalAtSite;
        };

        class GlobalLatticeData
        {
          public:
            GlobalLatticeData();
            ~GlobalLatticeData();

            void CollectFluidSiteDistribution();

            void SetBasicDetails(site_t iBlocksX,
                                 site_t iBlocksY,
                                 site_t iBlocksZ,
                                 site_t iBlockSize,
                                 distribn_t iVoxelSize,
                                 distribn_t iOriginX,
                                 distribn_t iOriginY,
                                 distribn_t iOriginZ);

            void GetThisRankSiteData(unsigned int *& bThisRankSiteData);

            site_t GetXSiteCount() const;
            site_t GetYSiteCount() const;
            site_t GetZSiteCount() const;
            site_t GetXBlockCount() const;
            site_t GetYBlockCount() const;
            site_t GetZBlockCount() const;

            site_t GetBlockSize() const;
            site_t GetBlockCount() const;

            distribn_t GetVoxelSize() const;
            distribn_t GetXOrigin() const;
            distribn_t GetYOrigin() const;
            distribn_t GetZOrigin() const;

            site_t GetSitesPerBlockVolumeUnit() const;

            bool IsValidBlock(site_t i, site_t j, site_t k) const;
            bool IsValidLatticeSite(site_t i, site_t j, site_t k) const;

            BlockData * Blocks;

            // Returns the type of collision/streaming update for the fluid site
            // with data "site_data".
            unsigned int GetCollisionType(unsigned int site_data) const;

            // Function that finds the pointer to the rank on which a particular site
            // resides. If the site is in an empty block, return NULL.
            const proc_t
            * GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const;

            // Function that gets the index of a block from its coordinates.
            site_t GetBlockIdFromBlockCoords(site_t blockI, site_t blockJ, site_t blockK) const;

            void GetBlockIJK(site_t block, site_t* i, site_t* j, site_t* k) const;
            site_t GetSiteCoord(site_t block, site_t localSiteCoord) const;
            const util::Vector3D<site_t>
            GetGlobalCoords(site_t blockNumber, const util::Vector3D<site_t>& localSiteCoords) const;
            unsigned int GetSiteData(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

            void ReadBlock(site_t block, io::writers::xdr::XdrReader* reader);
            const site_t* GetFluidSiteCountsOnEachProc() const;
            site_t GetFluidSiteCountOnProc(proc_t proc) const;

          public:
            // TODO public temporarily, until all usages are internal to the class.
            unsigned int log2BlockSize;
            // Array containing numbers of fluid sites on each processor.
            site_t* fluidSitesOnEachProcessor;
            // Hold the min and max site coordinates
            util::Vector3D<site_t> globalSiteMins, globalSiteMaxes;

          private:
            site_t mSitesPerBlockVolumeUnit;
            site_t mBlockCount;
            site_t mSitesX, mSitesY, mSitesZ;
            site_t mBlocksX, mBlocksY, mBlocksZ;
            site_t mBlockSize;
            distribn_t mVoxelSize;
            distribn_t mOriginX, mOriginY, mOriginZ;
        };

        class GeometryReader
        {
          public:
            GeometryReader(const bool reserveSteeringCore);
            ~GeometryReader();

            GlobalLatticeData* LoadAndDecompose(std::string& dataFilePath,
                                                reporting::Timers &timings);

          private:
            struct BlockLocation
            {
                site_t i, j, k;
            };

            void ReadPreamble(GlobalLatticeData* bGlobalLatticeData);

            void ReadHeader(site_t iBlockCount,
                            site_t* sitesInEachBlock,
                            unsigned int* bytesUsedByBlockInDataFile);

            void BlockDecomposition(const site_t iBlockCount,
                                    const GlobalLatticeData* iGlobLatDat,
                                    const site_t* fluidSitePerBlock,
                                    proc_t* initialProcForEachBlock);

            void DivideBlocks(site_t unassignedBlocks,
                              site_t totalBlockCount,
                              proc_t unitCount,
                              site_t* blocksOnEachUnit,
                              proc_t* unitForEachBlock,
                              const site_t* fluidSitesPerBlock,
                              const GlobalLatticeData* iGlobLatDat);

            void ReadInLocalBlocks(GlobalLatticeData* iGlobLatDat,
                                   const site_t* sitesPerBlock,
                                   const unsigned int* bytesPerBlock,
                                   const proc_t* unitForEachBlock,
                                   const proc_t localRank);

            void DecideWhichBlocksToRead(bool* readBlock,
                                         const proc_t* unitForEachBlock,
                                         const proc_t localRank,
                                         const GlobalLatticeData* iGlobLatDat);

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
            void ReadInBlock(GlobalLatticeData* iGlobLatDat,
                             MPI_Offset offsetSoFar,
                             char* buffer,
                             int* procsWantingThisBlockBuffer,
                             const site_t blockNumber,
                             const site_t sites,
                             const unsigned int bytes,
                             const int neededOnThisRank);

            /**
             * Calculates the number of the rank used to read in a given block.
             *
             * @param blockNumber
             * @return
             */
            proc_t GetReadingCoreForBlock(site_t blockNumber);

            bool Expand(std::vector<BlockLocation>* edgeBlocks,
                        std::vector<BlockLocation>* expansionBlocks,
                        const GlobalLatticeData* iGlobLatDat,
                        const site_t* fluidSitesPerBlock,
                        bool* blockAssigned,
                        proc_t currentUnit,
                        proc_t* unitForEachBlock,
                        site_t &blocksOnCurrentUnit,
                        site_t blocksPerUnit);

            void OptimiseDomainDecomposition(const site_t* sitesPerBlock,
                                             const unsigned int* bytesPerBlock,
                                             const proc_t* procForEachBlock,
                                             GlobalLatticeData* bGlobLatDat);

            void ValidateGraphData(idx_t* vtxDistribn,
                                   idx_t localVertexCount,
                                   idx_t* adjacenciesPerVertex,
                                   idx_t* adjacencies);

            void ValidateGlobLatDat(GlobalLatticeData* iGlobLatDat);

            void ValidateProcForEachBlock(proc_t* procForEachBlock, site_t blockCount);

            site_t GetHeaderLength(site_t blockCount) const;

            void GetSiteDistributionArray(idx_t* vertexDistribn,
                                          const site_t blockCount,
                                          const proc_t* procForEachBlock,
                                          const site_t* sitesPerBlock) const;

            void GetFirstSiteIndexOnEachBlock(idx_t* firstSiteIndexPerBlock,
                                              const site_t blockCount,
                                              const idx_t* vertexDistribution,
                                              const proc_t* procForEachBlock,
                                              const site_t* sitesPerBlock) const;

            void GetAdjacencyData(idx_t* adjacenciesPerVertex,
                                  idx_t* &localAdjacencies,
                                  const idx_t localVertexCount,
                                  const proc_t* procForEachBlock,
                                  const idx_t* firstSiteIndexPerBlock,
                                  const GlobalLatticeData* bGlobLatDat) const;

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
                                const idx_t* partitionVector,
                                const GlobalLatticeData* bGlobLatDat);

            void RereadBlocks(GlobalLatticeData* bGlobLatDat,
                              const idx_t* movesPerProc,
                              const idx_t* movesList,
                              const site_t* sitesPerBlock,
                              const unsigned int* bytesPerBlock,
                              const int* procForEachBlock);

            void ImplementMoves(GlobalLatticeData* bGlobLatDat,
                                const proc_t* procForEachBlock,
                                const idx_t* movesFromEachProc,
                                const idx_t* movesList) const;

            proc_t ConvertTopologyRankToGlobalRank(proc_t topologyRank) const;

            void CreateFileReadType(MPI_Datatype* dataType,
                                    const site_t blockCount,
                                    const bool* readBlock,
                                    const unsigned int* bytesPerBlock) const;

            // The config file starts with:
            // * 3 unsigned ints for the number of blocks in the x, y, z directions
            // * 1 unsigned int for the block size (number of sites along one edge of a block)
            // * 1 double for the voxel size
            // * 3 doubles for the world-position of site 0
            static const int preambleBytes = 4 * 4 + 4 * 8;
            static const proc_t HEADER_READING_RANK = 0;
            static const proc_t READING_GROUP_SIZE = 5;

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

        class NeighbouringProcessor
        {
          public:
            // Rank of the neighbouring processor.
            proc_t Rank;

            // The number of distributions shared between this neighbour and the current processor.
            site_t SharedFCount;

            // Index on this processor of the first distribution shared between this
            // neighbour and the current processor.
            site_t FirstSharedF;
        };

        /**
         * The protected default constructor does nothing. It exists to allow derivation from this
         * class for the purpose of testing.
         * @return
         */
        LatticeData();
        LatticeData(LocalLatticeData* localLattice, GlobalLatticeData* globalLattice);

        void SetSiteData(site_t siteIndex, unsigned int siteData);
        void SetWallNormal(site_t siteIndex, double normal[3]);
        void SetWallDistance(site_t siteIndex, double cutDistance[D3Q15::NUMVECTORS - 1]);
        void InitialiseNeighbourLookup(site_t** bSharedFLocationForEachProc,
                                       proc_t localRank,
                                       const unsigned int* iSiteDataForThisRank);
        void CountCollisionTypes(const unsigned int * lThisRankSiteData);
        void InitialisePointToPointComms(site_t** &lSharedFLocationForEachProc);
        void Initialise();

        LocalLatticeData localLatDat;
        GlobalLatticeData globLatDat;

        site_t totalFluidSites;

        site_t* f_recv_iv;
        // Number of local distributions shared with neighbouring processors.
        site_t totalSharedFs;
        // For each processor in the topology, holds the index into the
        // neighbouring processor vector.
        proc_t* neighbourIndexFromProcRank;
        // The vector of all neighbouring processors.
        std::vector<NeighbouringProcessor> neighbouringProcs;
    };
  }
}

#endif /* HEMELB_GEOMETRY_LATTICEDATA_H */
