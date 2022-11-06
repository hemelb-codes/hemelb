// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_DOMAIN_H
#define HEMELB_GEOMETRY_DOMAIN_H

#include <cassert>
#include <memory>
#include <vector>

#include "constants.h"
#include "units.h"
#include "geometry/Block.h"
#include "geometry/NeighbouringProcessor.h"
#include "geometry/Site.h"
#include "geometry/SiteDataBare.h"
#include "lb/lattices/LatticeInfo.h"
#include "reporting/Reportable.h"
#include "util/Vector3D.h"

namespace hemelb::extraction {
    class LocalDistributionInput;
}
namespace hemelb::lb {
    template<class TRAITS> class LBM; // IWYU pragma: keep
    struct InitialConditionBase;
}
namespace hemelb::net {
    class IOCommunicator;
}
namespace hemelb::reporting {
    class Dict;
}
namespace hemelb::tests::helpers {
    // Friend class to access all of domain_type's internals in tests
    class LatticeDataAccess;
}

namespace hemelb::geometry
{
    //class FieldData;
    class GmyReadResult;
    namespace neighbouring {
        class NeighbouringDomain;
    }

    // Hold geometrical and indexing type data about the domain to be simulated.
    class Domain : public reporting::Reportable
    {
        friend class FieldData;
        friend class extraction::LocalDistributionInput; //! Give access to the methods GetFOld and GetFNew.
        friend lb::InitialConditionBase;
        friend class tests::helpers::LatticeDataAccess;

      public:
        template<class TRAITS>
	friend class lb::LBM; //! Let the LBM have access to internals so it can initialise the distribution arrays.
        template<class LatticeData>
	friend class Site; //! Let the inner classes have access to site-related data that's otherwise private.

        Domain(const lb::lattices::LatticeInfo& latticeInfo, const GmyReadResult& readResult,
               const net::IOCommunicator& comms);

        ~Domain() override = default;

        inline net::IOCommunicator const & GetCommunicator() const
        {
          return comms;
        }

        /**
         * Get the lattice info object for the current lattice
         * @return
         */
        inline const lb::lattices::LatticeInfo& GetLatticeInfo() const
        {
          return latticeInfo;
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
        inline site_t GetLocalSiteIdFromLocalSiteCoords(
            const util::Vector3D<site_t>& siteCoords) const
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
        bool IsValidBlock(const util::Vector3D<site_t>& blockCoords) const;

        /**
         * Get the dimensions of the bounding box of the geometry in terms of blocks.
         * @return
         */
        inline const util::Vector3D<site_t>& GetBlockDimensions() const
        {
          return blockCounts;
        }

        /**
         * @param blockNumber
         * @return
         */
        inline const Block& GetBlock(site_t blockNumber) const
        {
          assert(blockNumber < site_t(blocks.size()));
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
        site_t GetContiguousSiteId(site_t x, site_t y, site_t z) const
        {
          return GetContiguousSiteId(LatticeVector(x, y, z));
        }

        /**
         * Get both the owner processor id and the local site id (if local)
         * @param globalLocation the location to retrieve information about
         * @param procId (out) the rank of the processor that owns the site
         * @param siteId (out) the index of the site for the property cache
         * @return true when globalLocation is local fluid, false otherwise
         */
        bool GetContiguousSiteId(const util::Vector3D<site_t>& globalLocation, proc_t& procId,
                                 site_t& siteId) const;

        /**
         * Get the global site coordinates from block coordinates and the site's local coordinates
         * within the block
         * @param blockCoords
         * @param localSiteCoords
         * @return
         */
        inline const util::Vector3D<site_t> GetGlobalCoords(
            const util::Vector3D<site_t>& blockCoords,
            const util::Vector3D<site_t>& localSiteCoords) const
        {
          return blockCoords * blockSize + localSiteCoords;
        }

        inline site_t GetGlobalNoncontiguousSiteIdFromGlobalCoords(
            const util::Vector3D<site_t>&globalCoords) const
        {
          return (globalCoords.x * sites.y + globalCoords.y) * sites.z + globalCoords.z;
        }

        inline site_t GetLocalContiguousIdFromGlobalNoncontiguousId(const site_t globalId) const
        {
          util::Vector3D<site_t> location;
          GetGlobalCoordsFromGlobalNoncontiguousSiteId(globalId, location);
          return GetContiguousSiteId(location);
        }

        void GetGlobalCoordsFromGlobalNoncontiguousSiteId(
            site_t globalId, util::Vector3D<site_t>& globalCoords) const
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

        void Report(reporting::Dict& dictionary) override;

        neighbouring::NeighbouringDomain &GetNeighbouringData();
        neighbouring::NeighbouringDomain const &GetNeighbouringData() const;

        int GetLocalRank() const;

        inline Site<Domain const> GetSite(site_t i) const {
            return Site<Domain const>{i, *this};
        }
      protected:
        /**
         * The protected default constructor does nothing. It exists to allow derivation from this
         * class for the purpose of testing.
         * @return
         */
        Domain(const lb::lattices::LatticeInfo& latticeInfo, const net::IOCommunicator& comms);

        void SetBasicDetails(util::Vector3D<site_t> blocks, site_t blockSize);

        void ProcessReadSites(const GmyReadResult& readResult);

        // TODO: should these parameters be copies?
        void PopulateWithReadData(
            const std::vector<site_t> midDomainBlockNumbers[COLLISION_TYPES],
            const std::vector<site_t> midDomainSiteNumbers[COLLISION_TYPES],
            const std::vector<SiteData> midDomainSiteData[COLLISION_TYPES],
            const std::vector<util::Vector3D<float> > midDomainWallNormals[COLLISION_TYPES],
            const std::vector<float> midDomainWallDistance[COLLISION_TYPES],
            const std::vector<site_t> domainEdgeBlockNumbers[COLLISION_TYPES],
            const std::vector<site_t> domainEdgeSiteNumbers[COLLISION_TYPES],
            const std::vector<SiteData> domainEdgeSiteData[COLLISION_TYPES],
            const std::vector<util::Vector3D<float> > domainEdgeWallNormals[COLLISION_TYPES],
            const std::vector<float> domainEdgeWallDistance[COLLISION_TYPES]);

        void CollectFluidSiteDistribution();
        void CollectGlobalSiteExtrema();

        void InitialiseNeighbourLookups();

        void InitialiseNeighbourLookup(
            std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialisePointToPointComms(
            std::vector<std::vector<site_t> >& sharedFLocationForEachProc);
        void InitialiseReceiveLookup(std::vector<std::vector<site_t> >& sharedFLocationForEachProc);

        sitedata_t GetSiteData(site_t iSiteI, site_t iSiteJ, site_t iSiteK) const;

        /**
         * Set the location (in the distribution array) of the streaming from site at siteIndex in direction
         * to be value.
         * @param iSiteIndex
         * @param iDirection
         * @param iValue
         */
        inline void SetNeighbourLocation(const site_t siteIndex, const unsigned int direction,
                                         const site_t distributionIndex)
        {
          neighbourIndices[siteIndex * latticeInfo.GetNumVectors() + direction] = distributionIndex;
        }

        void GetBlockIJK(site_t block, util::Vector3D<site_t>& blockCoords) const;

        // Method should remain protected, intent is to access this information via Site
        template<typename LatticeType>
        double GetCutDistance(site_t iSiteIndex, int iDirection) const
        {
          return distanceToWall[iSiteIndex * (LatticeType::NUMVECTORS - 1) + iDirection - 1];
        }

        /**
         * Get the wall normal at the given site
         * @param iSiteIndex
         * @return
         */
        // Method should remain protected, intent is to access this information via Site
        inline const util::Vector3D<distribn_t>& GetNormalToWall(site_t iSiteIndex) const
        {
          return wallNormalAtSite[iSiteIndex];
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

        /**
         * Get the site data object for the given index.
         * @param iSiteIndex
         * @return
         */
        // Method should remain protected, intent is to access this information via Site
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
        // Method should remain protected, intent is to access this information via Site
        inline SiteData &GetSiteData(site_t iSiteIndex)
        {
          return siteData[iSiteIndex];
        }

        // Method should remain protected, intent is to access this information via Site
        const distribn_t * GetCutDistances(site_t iSiteIndex) const
        {
          return &distanceToWall[iSiteIndex * (latticeInfo.GetNumVectors() - 1)];
        }

        distribn_t * GetCutDistances(site_t iSiteIndex)
        {
          return &distanceToWall[iSiteIndex * (latticeInfo.GetNumVectors() - 1)];
        }

        // Method should remain protected, intent is to access this information via Site
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
        util::Vector3D<site_t> sites;
        site_t sitesPerBlockVolumeUnit;
        site_t blockCount;

        site_t totalSharedFs; //! Number of local distributions shared with neighbouring processors.
        std::vector<NeighbouringProcessor> neighbouringProcs; //! Info about processors with neighbouring fluid sites.

        site_t midDomainProcCollisions[COLLISION_TYPES]; //! Number of fluid sites with all fluid neighbours on this rank, for each collision type.
        site_t domainEdgeProcCollisions[COLLISION_TYPES]; //! Number of fluid sites with at least one fluid neighbour on another rank, for each collision type.
        site_t localFluidSites; //! The number of local fluid sites.
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
        std::shared_ptr<neighbouring::NeighbouringDomain> neighbouringData;
        const net::IOCommunicator& comms;
    };

}


#endif /* HEMELB_GEOMETRY_DOMAIN_H */
