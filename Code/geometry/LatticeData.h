#ifndef HEMELB_GEOMETRY_LATTICEDATA_H
#define HEMELB_GEOMETRY_LATTICEDATA_H

#include <cstdio>
#include <vector>

#include "net/net.h"
#include "constants.h"
#include "configuration/SimConfig.h"
#include "geometry/Block.h"
#include "geometry/GeometryReader.h"
#include "geometry/Site.h"
#include "geometry/SiteData.h"
#include "reporting/Reportable.h"
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

    class LatticeData : public reporting::Reportable
    {
      public:
        friend class InnerSite<true> ;
        friend class InnerSite<false> ;
        friend class BlockTraverser;

        static LatticeData* Load(const bool reserveSteeringCore,
                                 const lb::lattices::LatticeInfo& latticeInfo,
                                 const std::string& dataFilePath,
                                 reporting::Timers &timings);

        virtual ~LatticeData();

        void SwapOldAndNew();
        void SendAndReceive(net::Net* net);
        void CopyReceived();

        const lb::lattices::LatticeInfo& GetLatticeInfo() const;

        Site GetSite(site_t localIndex);
        ConstSite GetSite(site_t localIndex) const;

        site_t GetXSiteCount() const;
        site_t GetYSiteCount() const;
        site_t GetZSiteCount() const;

        distribn_t GetVoxelSize() const;
        const util::Vector3D<distribn_t> GetOrigin() const;

        site_t GetBlockSize() const;
        site_t GetBlockCount() const;

        site_t GetSitesPerBlockVolumeUnit() const;

        /**
         * Convert a 3D site location within a block to a scalar site id on that block.
         *
         * @param blockCoords
         * @return
         */
        site_t GetLocalSiteIdFromLocalSiteCoords(const util::Vector3D<site_t>& siteCoords) const;
        site_t GetBlockIdFromBlockCoords(const util::Vector3D<site_t>& blockCoords) const;

        bool IsValidLatticeSite(const util::Vector3D<site_t>& siteCoords) const;

        distribn_t* GetFNew(site_t siteNumber);
        const distribn_t* GetFNew(site_t siteNumber) const;

        proc_t GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const;

        const Block& GetBlock(site_t blockNumber) const;

        site_t GetLocalFluidSiteCount() const;

        site_t GetContiguousSiteId(util::Vector3D<site_t> location) const;

        const util::Vector3D<site_t> GetGlobalCoords(const util::Vector3D<site_t>& blockCoords, const util::Vector3D<
            site_t>& localSiteCoords) const;

        const util::Vector3D<site_t>
        GetGlobalCoords(site_t blockNumber, const util::Vector3D<site_t>& localSiteCoords) const;
        util::Vector3D<site_t> GetSiteCoordsFromSiteId(site_t siteId) const;

        void GetBlockAndLocalSiteCoords(const util::Vector3D<site_t>& location,
                                        util::Vector3D<site_t>& blockCoords,
                                        util::Vector3D<site_t>& siteCoords) const;

        site_t GetMidDomainSiteCount() const;
        /**
         * Number of sites with all fluid neighbours residing on this rank, for the given
         * collision type.
         * @param collisionType
         * @return
         */
        site_t GetMidDomainCollisionCount(unsigned int collisionType) const;
        /**
         * Number of sites with at least one fluid neighbour residing on another rank
         * for the given collision type.
         * @param collisionType
         * @return
         */
        site_t GetDomainEdgeCollisionCount(unsigned int collisionType) const;

        site_t GetFluidSiteCountOnProc(proc_t proc) const;
        site_t GetTotalFluidSites() const;
        const util::Vector3D<site_t>& GetGlobalSiteMins() const;
        const util::Vector3D<site_t>& GetGlobalSiteMaxes() const;

        void Report(ctemplate::TemplateDictionary& dictionary);

      protected:
        /**
         * The protected default constructor does nothing. It exists to allow derivation from this
         * class for the purpose of testing.
         * @return
         */
        LatticeData(const lb::lattices::LatticeInfo& latticeInfo);
        LatticeData(const lb::lattices::LatticeInfo& latticeInfo, const Geometry& readResult);

        void SetBasicDetails(util::Vector3D<site_t> blocks, site_t blockSize, distribn_t voxelSize, util::Vector3D<
            distribn_t> originIn);

        void ProcessReadSites(const Geometry& readResult);

        void PopulateWithReadData(const std::vector<site_t> midDomainBlockNumbers[COLLISION_TYPES],
                                  const std::vector<site_t> midDomainSiteNumbers[COLLISION_TYPES],
                                  const std::vector<SiteData> midDomainSiteData[COLLISION_TYPES],
                                  const std::vector<util::Vector3D<float> > midDomainWallNormals[COLLISION_TYPES],
                                  const std::vector<float> midDomainWallDistance[COLLISION_TYPES],
                                  const std::vector<site_t> domainEdgeBlockNumbers[COLLISION_TYPES],
                                  const std::vector<site_t> domainEdgeSiteNumbers[COLLISION_TYPES],
                                  const std::vector<SiteData> domainEdgeSiteData[COLLISION_TYPES],
                                  const std::vector<util::Vector3D<float> > domainEdgeWallNormals[COLLISION_TYPES],
                                  const std::vector<float> domainEdgeWallDistance[COLLISION_TYPES])
        {
          // Populate the collision count arrays.
          for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
          {
            midDomainProcCollisions[collisionType] = midDomainBlockNumbers[collisionType].size();
            domainEdgeProcCollisions[collisionType] = domainEdgeBlockNumbers[collisionType].size();
          }
          // Data about local sites.
          localFluidSites = 0;
          // Data about contiguous local sites. First midDomain stuff, then domainEdge.
          for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
          {
            for (unsigned indexInType = 0; indexInType < midDomainProcCollisions[collisionType]; indexInType++)
            {
              siteData.push_back(midDomainSiteData[collisionType][indexInType]);
              wallNormalAtSite.push_back(midDomainWallNormals[collisionType][indexInType]);
              for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
              {
                distanceToWall.push_back(midDomainWallDistance[collisionType][indexInType
                    * (latticeInfo.GetNumVectors() - 1) + direction - 1]);
              }
              site_t blockId = midDomainBlockNumbers[collisionType][indexInType];
              site_t siteId = midDomainSiteNumbers[collisionType][indexInType];
              Blocks[blockId].SetLocalContiguousIndexForSite(siteId, localFluidSites);
              globalSiteCoords.push_back(GetGlobalCoords(blockId, GetSiteCoordsFromSiteId(siteId)));
              localFluidSites++;
            }

          }

          for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++)
          {
            for (unsigned indexInType = 0; indexInType < domainEdgeProcCollisions[collisionType]; indexInType++)
            {
              siteData.push_back(domainEdgeSiteData[collisionType][indexInType]);
              wallNormalAtSite.push_back(domainEdgeWallNormals[collisionType][indexInType]);
              for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
              {
                distanceToWall.push_back(domainEdgeWallDistance[collisionType][indexInType
                    * (latticeInfo.GetNumVectors() - 1) + direction - 1]);
              }
              site_t blockId = domainEdgeBlockNumbers[collisionType][indexInType];
              site_t siteId = domainEdgeSiteNumbers[collisionType][indexInType];
              Blocks[blockId].SetLocalContiguousIndexForSite(siteId, localFluidSites);
              globalSiteCoords.push_back(GetGlobalCoords(blockId, GetSiteCoordsFromSiteId(siteId)));
              localFluidSites++;
            }

          }

          fOld.resize(localFluidSites * latticeInfo.GetNumVectors() + 1 + totalSharedFs);
          fNew.resize(localFluidSites * latticeInfo.GetNumVectors() + 1 + totalSharedFs);
        }
        void CollectFluidSiteDistribution();
        void CollectGlobalSiteExtrema();

        void CleanEmptyBlocks();

        void InitialiseNeighbourStuff();

        void InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialisePointToPointComms(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        sitedata_t GetSiteData(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

        void SetNeighbourLocation(site_t iSiteIndex, unsigned int iDirection, site_t iValue);
        void GetBlockIJK(site_t block, util::Vector3D<site_t>& blockCoords) const;

        site_t GetXBlockCount() const;
        site_t GetYBlockCount() const;
        site_t GetZBlockCount() const;
        bool IsValidBlock(site_t i, site_t j, site_t k) const;

        const util::Vector3D<distribn_t>& GetNormalToWall(site_t iSiteIndex) const;
        distribn_t* GetFOld(site_t siteNumber);
        const distribn_t* GetFOld(site_t siteNumber) const;

        /*
         * This returns the index of the distribution to stream to.
         *
         * NOTE: If streaming would take the distribution out of the geometry, we instead stream
         * to the 'rubbish site', an extra position in the array that doesn't correspond to any
         * site in the geometry.
         */
        template<typename LatticeType>
        site_t GetStreamedIndex(site_t iSiteIndex, unsigned int iDirectionIndex) const
        {
          return neighbourIndices[iSiteIndex * LatticeType::NUMVECTORS + iDirectionIndex];
        }

        template<typename LatticeType>
        double GetCutDistance(site_t iSiteIndex, int iDirection) const
        {
          return distanceToWall[iSiteIndex * (LatticeType::NUMVECTORS - 1) + iDirection - 1];
        }

        SiteData GetSiteData(site_t iSiteIndex) const;

        /**
         * Get the global site coordinates from a contiguous site id.
         * @param siteIndex
         * @return
         */
        const util::Vector3D<site_t>& GetGlobalSiteCoords(site_t siteIndex) const;

        // Variables are listed here in approximate order of initialisation.
        // Note that all data is ordered in increasing order of collision type, by
        // midDomain (all neighbours on this core) then domainEdge (some neighbours on
        // another core)
        // I.e. midDomain type 0 to midDomain type 5 then domainEdge type 0 to domainEdge type 5.
        /**
         * Basic lattice variables.
         */
        const lb::lattices::LatticeInfo& latticeInfo;
        util::Vector3D<site_t> blockCounts;
        site_t blockSize;
        distribn_t voxelSize;
        util::Vector3D<distribn_t> origin;
        util::Vector3D<site_t> sites;
        site_t sitesPerBlockVolumeUnit;
        site_t blockCount;

        /**
         * Data about sending and receiving fs.
         */
        // Number of local distributions shared with neighbouring processors.
        site_t totalSharedFs;
        // Vector of info about processors with neighbouring fluid sites.
        std::vector<NeighbouringProcessor> neighbouringProcs;

        /**
         * Number of fluid sites with all fluid neighbours on this rank, for each collision type.
         */
        site_t midDomainProcCollisions[COLLISION_TYPES];
        /**
         * Number of fluid sites with at least one fluid neighbour on another rank, for each
         * collision type.
         */
        site_t domainEdgeProcCollisions[COLLISION_TYPES];

        site_t localFluidSites;

        std::vector<distribn_t> fOld;
        std::vector<distribn_t> fNew;

        /**
         * Data where local fluid sites are stored contiguously.
         */
        std::vector<Block> Blocks;

        std::vector<distribn_t> distanceToWall;
        /**
         * Hold the global site coordinates for each contiguous site.
         */
        std::vector<util::Vector3D<site_t> > globalSiteCoords;
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
