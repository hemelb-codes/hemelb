#ifndef HEMELB_GEOMETRY_LATTICEDATA_H
#define HEMELB_GEOMETRY_LATTICEDATA_H

#include <cstdio>
#include <vector>

#include "net/net.h"
#include "D3Q15.h"
#include "constants.h"
#include "configuration/SimConfig.h"
#include "mpiInclude.h"
#include "geometry/BlockTraverser.h"
#include "geometry/GeometryReader.h"
#include "geometry/SiteData.h"
// TODO Remove the stress type from the data file, so we can remove the dependence
// on LbmParams here.
#include "lb/LbmParameters.h"
#include "reporting/Timers.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
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

    class LatticeData
    {
      public:
        static LatticeData* Load(const bool reserveSteeringCore,
                                 std::string& dataFilePath,
                                 lb::LbmParameters* bLbmParams,
                                 reporting::Timers &timings);

        virtual ~LatticeData();

        void SwapOldAndNew();
        void SendAndReceive(net::Net* net);
        void CopyReceived();

        const util::Vector3D<distribn_t>& GetNormalToWall(site_t iSiteIndex) const;

        site_t GetXSiteCount() const;
        site_t GetYSiteCount() const;
        site_t GetZSiteCount() const;

        site_t GetXBlockCount() const;
        site_t GetYBlockCount() const;
        site_t GetZBlockCount() const;

        distribn_t GetVoxelSize() const;
        const util::Vector3D<distribn_t> GetOrigin() const;

        unsigned int GetLog2BlockSize() const;

        site_t GetBlockSize() const;
        site_t GetBlockCount() const;

        site_t GetSitesPerBlockVolumeUnit() const;

        site_t GetBlockIdFromBlockCoords(site_t i, site_t j, site_t k) const;

        bool IsValidBlock(site_t i, site_t j, site_t k) const;
        bool IsValidLatticeSite(site_t i, site_t j, site_t k) const;

        const proc_t
        * GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const;

        const BlockData* GetBlock(site_t blockNumber) const;
        BlockTraverser GetBlockTraverser() const;

        distribn_t* GetFOld(site_t siteNumber);
        distribn_t* GetFNew(site_t siteNumber);
        const distribn_t* GetFOld(site_t siteNumber) const;
        const distribn_t* GetFNew(site_t siteNumber) const;

        site_t GetLocalFluidSiteCount() const;
        SiteType GetSiteType(site_t iSiteIndex) const;
        int GetBoundaryId(site_t iSiteIndex) const;
        site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const;
        bool HasBoundary(const site_t iSiteIndex, const int iDirection) const;
        double GetCutDistance(site_t iSiteIndex, int iDirection) const;
        SiteData GetSiteData(site_t iSiteIndex) const;
        site_t GetContiguousSiteId(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;
        const util::Vector3D<site_t>
        GetGlobalCoords(site_t blockNumber, const util::Vector3D<site_t>& localSiteCoords) const;

        site_t GetInnerSiteCount() const;
        site_t GetInnerCollisionCount(unsigned int collisionType) const;
        site_t GetInterCollisionCount(unsigned int collisionType) const;
        unsigned int GetCollisionType(unsigned int site_data) const;
        const std::vector<site_t>& GetFluidSiteCountsOnEachProc() const;
        site_t GetFluidSiteCountOnProc(proc_t proc) const;
        site_t GetTotalFluidSites() const;
        const util::Vector3D<site_t>& GetGlobalSiteMins() const;
        const util::Vector3D<site_t>& GetGlobalSiteMaxes() const;

      protected:
        /**
         * The protected default constructor does nothing. It exists to allow derivation from this
         * class for the purpose of testing.
         * @return
         */
        LatticeData();
        LatticeData(const GeometryReadResult& readResult);

        void SetBasicDetails(util::Vector3D<site_t> blocks,
                             site_t blockSize,
                             distribn_t voxelSize,
                             util::Vector3D<distribn_t> originIn);

        void ProcessReadSites(const GeometryReadResult& readResult);

        void
            PopulateWithReadData(const std::vector<site_t> intraBlockNumbers[COLLISION_TYPES],
                                 const std::vector<site_t> intraSiteNumbers[COLLISION_TYPES],
                                 const std::vector<SiteData> intraSiteData[COLLISION_TYPES],
                                 const std::vector<util::Vector3D<double> > intraWallNormals[COLLISION_TYPES],
                                 const std::vector<double> intraWallDistance[COLLISION_TYPES],
                                 const std::vector<site_t> interBlockNumbers[COLLISION_TYPES],
                                 const std::vector<site_t> interSiteNumbers[COLLISION_TYPES],
                                 const std::vector<SiteData> interSiteData[COLLISION_TYPES],
                                 const std::vector<util::Vector3D<double> > interWallNormals[COLLISION_TYPES],
                                 const std::vector<double> interWallDistance[COLLISION_TYPES]);
        void CollectFluidSiteDistribution();
        void CollectGlobalSiteExtrema();

        void CleanEmptyBlocks();

        void InitialiseNeighbourStuff();

        void
        InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void
        InitialisePointToPointComms(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);

        void SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue);
        void GetBlockIJK(site_t block, site_t* blockI, site_t* blockJ, site_t* blockK) const;

        // Variables are listed here in approximate order of initialisation.
        // Note that all data is ordered in increasing order of collision type, by intra-proc then
        // inter-proc.
        // I.e. Intra type 0 to intra type 5 then inter type 0 to inter type 5.
        /**
         * Basic lattice variables.
         */
        util::Vector3D<site_t> blockCounts;
        site_t blockSize;
        distribn_t voxelSize;
        util::Vector3D<distribn_t> origin;
        util::Vector3D<site_t> sites;
        site_t sitesPerBlockVolumeUnit;
        // Log2 of the block size, useful for bit shifting.
        unsigned int log2BlockSize;
        site_t blockCount;

        /**
         * Data about sending and receiving fs.
         */
        // Number of local distributions shared with neighbouring processors.
        site_t totalSharedFs;
        // Vector of info about processors with neighbouring fluid sites.
        std::vector<NeighbouringProcessor> neighbouringProcs;

        /**
         * Data about local fluid sites.
         */
        site_t intraProcCollisions[COLLISION_TYPES];
        site_t interProcCollisions[COLLISION_TYPES];

        site_t localFluidSites;

        std::vector<distribn_t> fOld;
        std::vector<distribn_t> fNew;

        /**
         * Data where local fluid sites are stored contiguously.
         */
        std::vector<BlockData> Blocks;

        std::vector<distribn_t> distanceToWall;
        std::vector<util::Vector3D<distribn_t> > wallNormalAtSite;
        std::vector<SiteData> siteData;

        /**
         * Data about all fluid sites.
         */
        // Array containing numbers of fluid sites on each processor.
        std::vector<site_t> fluidSitesOnEachProcessor;
        site_t totalFluidSites;

        // Hold the min and max site coordinates
        util::Vector3D<site_t> globalSiteMins, globalSiteMaxes;

        /**
         * Data about neighbouring fluid sites.
         */
        std::vector<site_t> neighbourIndices;
        std::vector<site_t> f_recv_iv;

    };
  }
}

#endif /* HEMELB_GEOMETRY_LATTICEDATA_H */
