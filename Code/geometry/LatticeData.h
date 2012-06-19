#ifndef HEMELB_GEOMETRY_LATTICEDATA_H
#define HEMELB_GEOMETRY_LATTICEDATA_H

#include <cstdio>
#include <vector>

#include "net/net.h"
#include "constants.h"
#include "configuration/SimConfig.h"
#include "geometry/Block.h"
#include "geometry/GeometryReader.h"
#include "geometry/NeighbouringProcessor.h"
#include "geometry/Site.h"
#include "geometry/neighbouring/NeighbouringSite.h"
#include "geometry/SiteData.h"
#include "reporting/Reportable.h"
#include "reporting/Timers.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    class LatticeData : public reporting::Reportable
    {
      public:
        friend class BaseSite<LatticeData> ; //! Let the inner classes have access to site-related data that's otherwise private.
        friend class BaseSite<const LatticeData> ; //! Let the inner classes have access to site-related data that's otherwise private.

        LatticeData(const lb::lattices::LatticeInfo& latticeInfo, const Geometry& readResult);

        virtual ~LatticeData();

        /**
         * Swap the fOld and fNew arrays around.
         */
        inline void SwapOldAndNew()
        {
          oldDistributions.swap(newDistributions);
        }

        void SendAndReceive(net::Net* net);
        void CopyReceived();

        /**
         * Get the lattice info object for the current lattice
         * @return
         */
        inline const lb::lattices::LatticeInfo& GetLatticeInfo() const
        {
          return latticeInfo;
        }

        /**
         * Get a site object for the given index.
         * @param localIndex
         * @return
         */
        inline Site GetSite(site_t localIndex)
        {
          return Site(localIndex, *this);
        }

        /**
         * Get a const site object for the given index.
         * @param localIndex
         * @return
         */
        inline ConstSite GetSite(site_t localIndex) const
        {
          return ConstSite(localIndex, *this);
        }

        /**
         * Get the dimensions of the bounding-box of the geometry in sites
         * @return
         */
        inline const util::Vector3D<site_t>& GetSiteDimensions() const
        {
          return sites;
        }

        /**
         * Get the voxel size, the real-world length of a lattice unit
         * @return
         */
        inline distribn_t GetVoxelSize() const
        {
          return voxelSize;
        }

        /**
         * Get the real-world origin of the geometry (i.e. where the lattice point
         * (0,0,0) maps to)
         * @return
         */
        inline const util::Vector3D<distribn_t> GetOrigin() const
        {
          return origin;
        }

        /**
         * Get the number of sites along one block length
         * @return
         */
        inline site_t GetBlockSize() const
        {
          return blockSize;
        }

        /**
         * Get the total number of blocks in the bounding-box.
         * @return
         */
        inline site_t GetBlockCount() const
        {
          return blockCount;
        }

        /**
         * Get the total number of sites in a block-cube.
         * @return
         */
        inline site_t GetSitesPerBlockVolumeUnit() const
        {
          return sitesPerBlockVolumeUnit;
        }

        /**
         * Convert a 3D site location within a block to a scalar site id on that block.
         *
         * @param blockCoords
         * @return
         */
        inline site_t GetLocalSiteIdFromLocalSiteCoords(const util::Vector3D<site_t>& siteCoords) const
        {
          return ( (siteCoords.x * blockSize) + siteCoords.y) * blockSize + siteCoords.z;
        }

        /**
         * Get the scalar id of a block from its 3D coordinates
         * @param blockCoords
         * @return
         */
        inline site_t GetBlockIdFromBlockCoords(const util::Vector3D<site_t>& blockCoords) const
        {
          return (blockCoords.x * blockCounts.y + blockCoords.y) * blockCounts.z + blockCoords.z;
        }

        bool IsValidLatticeSite(const util::Vector3D<site_t>& siteCoords) const;

        /**
         * Get a pointer into the fNew array at the given index
         * @param distributionIndex
         * @return
         */
        inline distribn_t* GetFNew(site_t distributionIndex)
        {
          return &newDistributions[distributionIndex];
        }

        /**
         * Get a pointer into the fNew array at the given index. This version of the above lets us
         * use a const version of a LatticeData to get a const *.
         * @param distributionIndex
         * @return
         */
        inline const distribn_t* GetFNew(site_t siteNumber) const
        {
          return &newDistributions[siteNumber];
        }

        proc_t GetProcIdFromGlobalCoords(const util::Vector3D<site_t>& globalSiteCoords) const;

        /**
         * True if the given coordinates correspond to a valid block within the bounding
         * box of the geometry.
         *
         * @param i
         * @param j
         * @param k
         * @return
         */
        bool IsValidBlock(site_t i, site_t j, site_t k) const;

        /**
         * Get the dimensions of the bounding box of the geometry in terms of blocks.
         * @return
         */
        inline const util::Vector3D<site_t>& GetBlockDimensions() const
        {
          return blockCounts;
        }

        /**
         * Get the block data for the given block id.
         * @param blockNumber
         * @return
         */
        inline const Block& GetBlock(site_t blockNumber) const
        {
          return blocks[blockNumber];
        }

        /**
         * Get the number of fluid sites local to this proc.
         * @return
         */
        inline site_t GetLocalFluidSiteCount() const
        {
          return localFluidSites;
        }

        site_t GetContiguousSiteId(util::Vector3D<site_t> location) const;

        /**
         * Get the global site coordinates from block coordinates and the site's local coordinates
         * within the block
         * @param blockCoords
         * @param localSiteCoords
         * @return
         */
        inline const util::Vector3D<site_t> GetGlobalCoords(const util::Vector3D<site_t>& blockCoords,
                                                            const util::Vector3D<site_t>& localSiteCoords) const
        {
          return blockCoords * blockSize + localSiteCoords;
        }

        inline site_t GetGlobalNoncontiguousSiteIdFromGlobalCoords(const util::Vector3D<site_t>&globalCoords) const
        {
          return (globalCoords.x * sites.y + globalCoords.y) * sites.z + globalCoords.z;
        }

        inline site_t GetLocalContiguousIdFromGlobalNoncontiguousId(const site_t globalId) const
        {
          util::Vector3D<site_t> location;
          GetGlobalCoordsFromGlobalNoncontiguousSiteId(globalId, location);
          return GetContiguousSiteId(location);
        }

        void GetGlobalCoordsFromGlobalNoncontiguousSiteId(site_t globalId, util::Vector3D<site_t>& globalCoords) const
        {
          globalCoords.z = globalId % sites.z;
          site_t blockIJData = globalId / sites.z;
          globalCoords.y = blockIJData % sites.y;
          globalCoords.x = blockIJData / sites.y;
        }

        proc_t ProcProvidingSiteByGlobalNoncontiguousId(site_t globalId) const
        {
          util::Vector3D<site_t> resultCoord;
          GetGlobalCoordsFromGlobalNoncontiguousSiteId(globalId, resultCoord);
          return GetProcIdFromGlobalCoords(resultCoord);
        }

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
        inline site_t GetMidDomainCollisionCount(unsigned int collisionType) const
        {
          return midDomainProcCollisions[collisionType];
        }

        /**
         * Number of sites with at least one fluid neighbour residing on another rank
         * for the given collision type.
         * @param collisionType
         * @return
         */
        inline site_t GetDomainEdgeCollisionCount(unsigned int collisionType) const
        {
          return domainEdgeProcCollisions[collisionType];
        }

        /**
         * Get the number of fluid sites on the given rank.
         * @param proc
         * @return
         */
        inline site_t GetFluidSiteCountOnProc(proc_t proc) const
        {
          return fluidSitesOnEachProcessor[proc];
        }

        /**
         * Get the total number of fluid sites in the whole geometry.
         * @return
         */
        inline site_t GetTotalFluidSites() const
        {
          return totalFluidSites;
        }

        /**
         * Get the minimal x,y,z coordinates for any fluid site.
         * @return
         */
        inline const util::Vector3D<site_t>& GetGlobalSiteMins() const
        {
          return globalSiteMins;
        }

        /**
         * Get the maximal x,y,z coordinates for any fluid site.
         * @return
         */
        inline const util::Vector3D<site_t>& GetGlobalSiteMaxes() const
        {
          return globalSiteMaxes;
        }

        void Report(ctemplate::TemplateDictionary& dictionary);

        neighbouring::NeighbouringLatticeData &GetNeighbouringData();
        neighbouring::NeighbouringLatticeData const &GetNeighbouringData() const;
      protected:
        /**
         * The protected default constructor does nothing. It exists to allow derivation from this
         * class for the purpose of testing.
         * @return
         */
        LatticeData(const lb::lattices::LatticeInfo& latticeInfo);

        void SetBasicDetails(util::Vector3D<site_t> blocks,
                             site_t blockSize,
                             distribn_t voxelSize,
                             util::Vector3D<distribn_t> originIn);

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
              blocks[blockId].SetLocalContiguousIndexForSite(siteId, localFluidSites);
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
              blocks[blockId].SetLocalContiguousIndexForSite(siteId, localFluidSites);
              globalSiteCoords.push_back(GetGlobalCoords(blockId, GetSiteCoordsFromSiteId(siteId)));
              localFluidSites++;
            }

          }

          oldDistributions.resize(localFluidSites * latticeInfo.GetNumVectors() + 1 + totalSharedFs);
          newDistributions.resize(localFluidSites * latticeInfo.GetNumVectors() + 1 + totalSharedFs);
        }
        void CollectFluidSiteDistribution();
        void CollectGlobalSiteExtrema();

        void CleanEmptyBlocks();

        void InitialiseNeighbourLookups();

        void InitialiseNeighbourLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialisePointToPointComms(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);

        sitedata_t GetSiteData(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

        /**
         * Set the location (in the distribution array) of the streaming from site at siteIndex in direction
         * to be value.
         * @param iSiteIndex
         * @param iDirection
         * @param iValue
         */
        inline void SetNeighbourLocation(const site_t siteIndex,
                                         const unsigned int direction,
                                         const site_t distributionIndex)
        {
          neighbourIndices[siteIndex * latticeInfo.GetNumVectors() + direction] = distributionIndex;
        }

        void GetBlockIJK(site_t block, util::Vector3D<site_t>& blockCoords) const;

        /**
         * Get the wall normal at the given site
         * @param iSiteIndex
         * @return
         */
        inline const util::Vector3D<distribn_t>& GetNormalToWall(site_t iSiteIndex) const
        {
          return wallNormalAtSite[iSiteIndex];
        }

        /**
         * Get a pointer to the fOld array starting at the requested index
         * @param distributionIndex
         * @return
         */
        distribn_t* GetFOld(site_t distributionIndex)
        {
          return &oldDistributions[distributionIndex];
        }

        /**
         * Get a pointer to the fOld array starting at the requested index. This version
         * of the function allows us to access the fOld array in a const way from a const
         * LatticeData
         * @param distributionIndex
         * @return
         */
        const distribn_t* GetFOld(site_t distributionIndex) const
        {
          return &oldDistributions[distributionIndex];
        }

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

        /**
         * Get the site data object for the given index.
         * @param iSiteIndex
         * @return
         */
        inline const SiteData &GetSiteData(site_t iSiteIndex) const
        {
          return siteData[iSiteIndex];
        }

        /***
         * Non-const version of getting site-data, for use with MPI calls, where const-ness is not respected on sends.
         * Not available on const LatticeDatas
         * @param iSiteIndex
         * @return
         */
        inline SiteData &GetSiteData(site_t iSiteIndex)
        {
          return siteData[iSiteIndex];
        }

        distribn_t * GetCutDistances(site_t iSiteIndex)
        {
          return &distanceToWall[iSiteIndex * (latticeInfo.GetNumVectors() - 1)];
        }

        util::Vector3D<distribn_t>& GetNormalToWall(site_t iSiteIndex)
        {
          return wallNormalAtSite[iSiteIndex];
        }
        /**
         * Get the global site coordinates from a contiguous site id.
         * @param siteIndex
         * @return
         */
        inline const util::Vector3D<site_t>& GetGlobalSiteCoords(site_t siteIndex) const
        {
          return globalSiteCoords[siteIndex];
        }

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

        site_t totalSharedFs; //! Number of local distributions shared with neighbouring processors.
        std::vector<NeighbouringProcessor> neighbouringProcs; //! Info about processors with neighbouring fluid sites.

        site_t midDomainProcCollisions[COLLISION_TYPES]; //! Number of fluid sites with all fluid neighbours on this rank, for each collision type.
        site_t domainEdgeProcCollisions[COLLISION_TYPES]; //! Number of fluid sites with at least one fluid neighbour on another rank, for each collision type.
        site_t localFluidSites; //! The number of local fluid sites.
        std::vector<distribn_t> oldDistributions; //! The distribution values for the previous time step.
        std::vector<distribn_t> newDistributions; //! The distribution values for the next time step.
        std::vector<Block> blocks; //! Data where local fluid sites are stored contiguously.

        std::vector<distribn_t> distanceToWall; //! Hold the distance to the wall for each fluid site.
        std::vector<util::Vector3D<site_t> > globalSiteCoords; //! Hold the global site coordinates for each contiguous site.
        std::vector<util::Vector3D<distribn_t> > wallNormalAtSite; //! Holds the wall normal near the fluid site, where appropriate
        std::vector<SiteData> siteData; //! Holds the SiteData for each site.
        std::vector<site_t> fluidSitesOnEachProcessor; //! Array containing numbers of fluid sites on each processor.
        site_t totalFluidSites; //! The total number of fluid sites in the geometry.
        util::Vector3D<site_t> globalSiteMins, globalSiteMaxes; //! The minimal and maximal coordinates of any fluid sites.
        std::vector<site_t> neighbourIndices; //! Data about neighbouring fluid sites.
        std::vector<site_t> streamingIndicesForReceivedDistributions; //! The indices to stream to for distributions received from other processors.
        neighbouring::NeighbouringLatticeData *neighbouringData;
    };
  }
}

#endif /* HEMELB_GEOMETRY_LATTICEDATA_H */
